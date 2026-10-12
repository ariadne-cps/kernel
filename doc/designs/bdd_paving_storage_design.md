# BDD-backed paving storage design

This document is a working design for introducing an ROBDD-backed grid-paving storage implementation in `ariadne-cps/kernel` while retaining the existing explicit binary-tree implementation.

The goal is not to replace `GridTreePaving` immediately. The two storage representations should coexist behind the paving abstraction so that they can be compared for correctness, memory use, sequential performance and concurrent scalability on the same geometric workloads.

The decision-diagram implementation belongs in `ariadne-cps/bdd`. Kernel owns the geometric interpretation that maps decision levels to grid-cell subdivision paths.

## 1. Goals

The paving work has five primary goals.

1. Preserve the existing tree-backed paving implementation as a supported representation and performance baseline.
2. Add an ROBDD-backed paving implementation with equivalent geometric semantics.
3. Make the existing paving interfaces genuinely representation-independent so algorithms can consume either backend.
4. Use the tree implementation as a differential correctness oracle for the ROBDD implementation.
5. Drive the missing ROBDD structural API from the requirements of real paving operations, including the historical BuDDy-backed implementation in `ariadne-cps/release-1.0`.

A secondary goal is to make performance comparisons meaningful. The same high-level operation should be exercisable through `PavingInterface` against both representations without changing the geometric algorithm.

## 2. Existing Kernel abstraction

Kernel already contains the correct high-level abstraction boundary.

`SubPavingInterface` provides the non-owning/logical subpaving view and query surface. It currently includes:

- `dimension()`;
- `grid()`;
- `size()`;
- `root_cell()`;
- `superset(GridCell)`;
- `equals(SubPavingInterface)`;
- `subset(SubPavingInterface)`;
- `intersects(SubPavingInterface)`;
- `branch(Bool)`;
- forward iteration through enabled `GridCell` values;
- `bounding_box()`;
- `inside()`, `separated()` and `overlaps()`;
- `set_root_cell()`;
- `mince()`;
- `recombine()`.

`PavingInterface` extends `SubPavingInterface` with mutating paving operations, including:

- `smallest_enclosing_primary_cell()`;
- outer and inner approximation insertion;
- `adjoin(SubPavingInterface)`;
- `restrict(SubPavingInterface)`;
- `remove(SubPavingInterface)`;
- cell insertion/removal.

The intended representation hierarchy is therefore:

```text
SubPavingInterface
    +-- GridTreeSubpaving
    +-- GridRobddSubpaving

PavingInterface
    +-- GridTreePaving
    +-- GridRobddPaving
```

A new storage-specific abstract interface is not currently justified. The existing paving interfaces should instead be made truly backend-independent.

## 3. Current representation leak

The existing interface is only partially representation-independent.

Several `GridTreeSubpaving` operations accept `SubPavingInterface const&` but immediately require the operand to be a `GridTreeSubpaving`:

```cpp
Bool GridTreeSubpaving::equals(const SubPavingInterface& paving) const {
    return (*this) == dynamic_cast<const GridTreeSubpaving&>(paving);
}

Bool GridTreeSubpaving::subset(const SubPavingInterface& paving) const {
    return Ariadne::subset(
        *this,
        dynamic_cast<const GridTreeSubpaving&>(paving));
}

Bool GridTreeSubpaving::intersects(const SubPavingInterface& paving) const {
    return Ariadne::intersect(
        *this,
        dynamic_cast<const GridTreeSubpaving&>(paving));
}
```

The same issue appears in the mutating operations:

```cpp
GridTreePaving::adjoin(SubPavingInterface const&)
GridTreePaving::restrict(SubPavingInterface const&)
GridTreePaving::remove(SubPavingInterface const&)
```

These currently perform a mandatory cast to `GridTreeSubpaving`. A mixed operation such as

```text
GridTreePaving op GridRobddPaving
```

would therefore fail rather than use the common paving semantics.

This is the first abstraction issue to resolve.

## 4. Representation-independent operations

The implementation should distinguish homogeneous fast paths from heterogeneous fallbacks.

Conceptually:

```text
same representation?
    |
    +-- yes -> representation-specific fast path
    |
    +-- no  -> generic PavingInterface/SubPavingInterface fallback
```

For example, tree/tree insertion can retain its existing recursive tree merge, while ROBDD/ROBDD insertion can be a Boolean disjunction of canonical roots.

A mixed operation can initially use the public cell semantics:

```cpp
Void GridTreePaving::adjoin(SubPavingInterface const& other) {
    if (auto tree = dynamic_cast<GridTreeSubpaving const*>(&other)) {
        this->adjoin(*tree);
        return;
    }

    for (GridCell const& cell : other) {
        this->adjoin(cell);
    }
}
```

The exact implementation need not have this shape, but the cast must become an optimization choice rather than a precondition.

For `equals`, `subset` and `intersects`, a generic fallback needs more care because blindly enumerating all cells may be expensive. Correctness comes before optimization during initial integration. Once both representations exist, cross-representation algorithms can be profiled and specialized only where measurements justify them.

Double dispatch is not required initially. With two concrete representations, optional same-type fast paths plus generic fallback keep the design simpler. A richer dispatch mechanism can be introduced later if the number of paving backends grows.

## 5. Storage ownership boundary

The dependency direction must remain:

```text
ariadne-cps/bdd
      ^
      |
Kernel geometry
```

The BDD repository must not know about:

- `Grid`;
- `GridCell`;
- `PavingInterface`;
- root extents;
- geometric dimension;
- Ariadne set interfaces.

Kernel should provide the geometric adapter that interprets ROBDD variable levels as subdivision levels.

The two concrete Kernel storage paths are intended to look like:

```text
GridTreePaving
    |
    +-- BinaryTreeNode

GridRobddPaving
    |
    +-- Robdd root
    +-- shared RobddManager/context
```

`BinaryTreeNode` remains a Kernel geometry implementation detail. It should not be moved into the BDD repository merely because it is an ordered binary decision tree.

## 6. Proposed ROBDD paving representation

A first `GridRobddPaving` is expected to contain conceptually:

```text
Grid
primary/root extent
logical mince depth
ROBDD context/manager
ROBDD root
```

A `GridRobddSubpaving` is a view into the same logical paving. It needs enough information to interpret the current root cell and the corresponding logical decision level without materializing a physical tree node.

An indicative model is:

```text
GridRobddSubpaving
    grid
    root extent
    path/root-cell metadata
    logical next decision level
    immutable Robdd root
    shared manager/context
```

The exact fields should be derived from the existing `GridTreeSubpaving` semantics and the historical `BDDTreeSet`, not fixed prematurely by this document.

### 6.1 Skipped levels

A reduced ROBDD may omit subdivision variables.

For grid storage, an omitted level means that both children have identical represented content. Therefore a logical paving branch may exist even when there is no physical ROBDD node at that level.

Conceptually:

```text
logical level L
    |
    +-- ROBDD top variable == L
    |       branch(0) -> low
    |       branch(1) -> high
    |
    +-- ROBDD top variable > L or terminal
            branch(0) -> same root
            branch(1) -> same root
```

Traversal, branching, enabled-cell counting and iteration must all use logical levels rather than assuming a physical ROBDD node exists for every subdivision.

### 6.2 Terminals

`false` represents a uniformly disabled logical subtree.

`true` represents a uniformly enabled logical subtree through the currently required paving depth. Iteration and `size()` must account for the number of implicit cells represented by a terminal root when the logical/mince depth extends below it.

This is a major semantic difference from simply counting ROBDD nodes.

## 7. Manager/context ownership

Pavings that participate in efficient Boolean set algebra must share the same canonical ROBDD universe.

A manager per paving would make

```text
A union B
A intersection B
A difference B
```

require cross-manager import/reconstruction before every operation.

The preferred architecture is:

```text
GridRobddContext
        |
        +-- shared fixed-order RobddManager
        |
        +-- paving A -> immutable root A
        +-- paving B -> immutable root B
        +-- paving C -> immutable root C
```

The exact scope of the context must be decided during implementation. It may be explicit in the paving API or supplied by an internal/shared service, but manager compatibility must never be accidental.

ROBDD nodes remain immutable. A mutable paving operation computes a new root and replaces the paving's root handle. Copies can therefore share the manager and canonical node storage.

This creates an important performance contrast with `GridTreePaving`, whose current copy semantics recursively duplicate the explicit tree. Copy and snapshot workloads should be included in profiling.

## 8. Variable ordering for paving

Grid paving requires fixed natural variable ordering because the variable sequence is the subdivision path semantics.

Dynamic ROBDD variable reordering must therefore be unavailable for the manager/context used by `GridRobddPaving`.

The BDD library may support reorderable ROBDD managers for SAT/general-purpose use, but the paving backend must select the fixed-order specialization and must not pay repeated runtime checks for a feature it cannot use.

Logical BDD variable identity and ordering level should remain distinct concepts in the general BDD design. In the paving fixed-order specialization they may coincide as an optimized representation.

## 9. Mince and recombine semantics

The existing tree and ROBDD representations implement these concepts differently.

For the explicit tree:

```text
mince
    -> physically split leaves to the requested depth

recombine
    -> physically remove equal sibling leaves
```

For a reduced ROBDD:

```text
mince
    -> increase logical interpretation depth
       without necessarily allocating nodes

recombine
    -> reduction is already canonical;
       discard artificial logical mince depth where appropriate
```

The historical BuDDy-backed paving followed the latter pattern.

Therefore the common interface must specify `mince` and `recombine` semantically rather than in terms of node creation/deletion. Even if these methods are later deprecated, the migration phase needs equivalent observable behavior in both backends so that historical tests remain meaningful.

## 10. Branch and subpaving semantics

`SubPavingInterface::branch(Bool)` already provides a suitable storage-independent operation.

For the tree backend, a branch is naturally represented by a child `BinaryTreeNode*`.

For the ROBDD backend, a branch is a logical cofactor/view at the next subdivision level. If that level is absent in the reduced diagram, both branches retain the same ROBDD root while advancing the logical path metadata.

This property is useful because it allows existing algorithms written recursively in terms of subpavings to remain independent of whether a concrete node exists at every level.

## 11. Same-backend fast paths

The initial optimized mapping should be:

| Paving operation | Tree/Tree | ROBDD/ROBDD |
| --- | --- | --- |
| adjoin / union | recursive tree merge | `OR` |
| restrict / intersection | recursive tree restriction | `AND` |
| remove / difference | recursive tree removal | `A AND NOT B` |
| equality | structural/semantic tree comparison | canonical root equality after compatible alignment |
| subset | tree recursion | Boolean implication test, e.g. `A AND NOT B == false` or specialized primitive |
| intersects | tree recursion | `A AND B != false` or specialized primitive |
| branch | child node | logical cofactor/view |
| recombine | physical sibling collapse | canonical reduction / logical-depth adjustment |

This table describes semantic direction, not a requirement that every operation be implemented by materializing an intermediate ROBDD. Specialized predicates such as disjointness or subset may later avoid constructing a full result if profiling shows a benefit.

## 12. Heterogeneous operations

Tree/ROBDD and ROBDD/tree operations must be correct even when they are not initially optimal.

The simplest correctness fallback is conversion through enabled grid cells when the root grids and cell semantics are compatible.

Potential future strategies include:

- bulk import of tree paths into an ROBDD;
- ROBDD traversal directly against a tree without intermediate conversion;
- temporary conversion of the smaller operand;
- representation-aware double dispatch.

No cross-backend optimization should be added before profiling demonstrates that the generic fallback materially affects real workloads.

Cross-backend tests are mandatory because they verify that `PavingInterface` is actually an abstraction rather than a common base class that only works within one concrete type.

## 13. Existing Kernel coupling to the tree backend

In addition to mandatory casts inside `GridTreePaving`, Kernel currently contains concrete-type assumptions such as:

```cpp
typedef GridTreePaving SetApproximationType;
```

in set-interface-related code.

These assumptions should be catalogued before switching any default representation.

The first milestone does not require changing the public/default paving type. It is preferable to:

1. make `PavingInterface` genuinely representation-independent;
2. implement `GridRobddPaving`;
3. establish differential correctness and performance data;
4. only then decide whether the default approximation type, factory mechanism or user-facing aliases should change.

Possible later selection mechanisms include an explicit factory or storage policy, but the choice should be driven by actual usage rather than introduced speculatively.

## 14. ROBDD extensions required by paving

The current ROBDD API supports Boolean construction and evaluation but is not yet sufficient for the historical grid-paving algorithms.

The historical comparison between the explicit tree representation and the BuDDy-backed `BDDTreeSet` in `ariadne-cps/release-1.0` identified the following required capabilities.

### 14.1 Structural inspection

Kernel needs to inspect the canonical diagram without exposing the BDD arena itself.

Required operations include:

- terminal test;
- top variable of a non-terminal diagram;
- low cofactor;
- high cofactor;
- possibly a combined node-inspection operation returning top/low/high together to avoid repeated arena lookups.

A combined inspection API is especially relevant because the current ROBDD implementation serializes arena reads through a mutex. Paving traversal should not acquire a lock independently for `top`, `low` and `high` if one coherent read can provide all three.

The public API should preserve manager provenance checks at its boundary and then perform recursive traversal using internal node references/identifiers without repeated public-handle validation.

### 14.2 Construction from a decision variable and branches

Historical tree/BDD transformations need the ability to construct a canonical decision node from:

```text
(variable, low, high)
```

while preserving the ROBDD reduction and unique-table invariants.

This must not expose raw unique-table mutation. A public structural constructor or an equivalent higher-level primitive should delegate to the manager's canonical node factory.

### 14.3 Support extraction

Projection and variable transformations need the set of variables actually present in a diagram.

The new ROBDD layer therefore needs support extraction analogous in purpose to BuDDy's support/scanset operations.

The result should use Ariadne BDD variable identities, not expose implementation ordering levels as if they were stable logical identities.

### 14.4 Variable shift and renaming/replacement

Historical paving algorithms change root height and dimension mappings by shifting or replacing BDD variables.

Required functionality includes:

- efficient order-preserving shift of a range of variables;
- general variable renaming/replacement where the resulting ordering is valid;
- explicit validation of mappings that would violate the manager's ordering invariant.

Root-height equalization, height increase/minimization and projection all depend on this functionality.

The API should express Ariadne semantics directly rather than mimic BuDDy's pair-table API.

### 14.5 Existential quantification

Projection removes dimensions/decision variables. The ROBDD manager therefore needs existential quantification over a selected variable set.

Conceptually:

```text
exists(vars, f)
```

must return the canonical ROBDD representing existential elimination of those variables.

The implementation should support the fixed-order paving manager efficiently and remain usable by later SAT/model-set applications.

### 14.6 Fused conjunction plus existential quantification

The historical projection path uses BuDDy's:

```text
bdd_appex(lhs, rhs, bddop_and, vars)
```

which combines conjunction with existential quantification.

Ariadne BDD should provide an optimized equivalent, tentatively:

```text
and_exists(lhs, rhs, vars)
```

rather than requiring:

```text
exists(vars, conjunction(lhs, rhs))
```

to materialize a potentially much larger intermediate diagram.

This operation is a high-priority paving primitive because multidimensional projection can otherwise create avoidable memory and cache pressure.

### 14.7 Efficient traversal with skipped variables

The ROBDD library should expose enough structural information for Kernel to distinguish:

```text
top variable == requested logical level
top variable > requested logical level
terminal
```

without expanding the reduced diagram.

This is required for:

- branch views;
- enabled-cell iteration;
- cell counting;
- subset/disjoint traversal;
- conversion to explicit tree paths;
- debugging and differential comparison.

The library should not itself know the paving's logical maximum depth; Kernel owns interpretation of missing levels.

### 14.8 Fixed-order manager capability

The paving manager must provide a configuration in which dynamic reordering is structurally unavailable.

The expected initial paving configuration is:

```text
ROBDD
fixed natural ordering
regular edges
```

Complemented edges are not required for the first grid-paving implementation.

They should later be benchmarked because difference/remove may benefit from cheap complement and increased sharing, while traversal-heavy operations may pay a small polarity-propagation cost.

The fixed + regular-edge path should not contain repeated hot-path feature checks for optional reordering or complemented edges.

### 14.9 Node identity versus node reference

Before complemented edges are introduced, the BDD implementation should distinguish:

```text
NodeId
    physical canonical node identity

NodeRef
    semantic reference to a node, potentially carrying edge metadata
```

This separation is primarily an internal BDD concern but matters to Kernel because structural traversal must operate on semantic ROBDD values, not accidentally interpret a physical node without applying reference polarity.

### 14.10 Shared-manager scalability

Paving is expected to keep many roots alive in one manager/context. Therefore scalability of one shared manager is more important than optimizing many isolated managers.

The current single mutex and ordered maps are a correctness baseline only.

Before claiming a concurrent paving advantage, BDD needs:

- concurrent-readable immutable node storage;
- fewer/reduced locks on structural reads;
- sharded hash-based unique table;
- bounded/lossy hash-based computed table where appropriate;
- instrumentation for table hits/misses and node allocation;
- structured parallel recursion using `ariadne-cps/threading`;
- ThreadSanitizer stress tests covering BDD-owned shared state.

### 14.11 Lifetime and garbage collection

Long-running paving algorithms may create many temporary roots while retaining only a smaller live set.

The append-only manager is sufficient for early correctness/performance comparison, but memory-growth measurements must be recorded from the beginning.

Garbage collection or another reclamation strategy should be introduced only after node publication and concurrent lifetime rules are stable. Premature concurrent reclamation would complicate the first paving integration substantially.

## 15. Historical release-1.0 behavior to preserve

The migration target is semantic compatibility with the historical tree and BuDDy implementations, not source compatibility with the BuDDy API.

The historical behavior that must be represented in compatibility tests includes:

- empty and full paving states;
- root-height increase;
- root-height minimization;
- root-cell alignment/equalization before binary operations;
- mince/recombine behavior;
- union/adjoin;
- intersection/restrict;
- difference/remove;
- subset;
- disjoint/intersection tests;
- enabled-cell counting;
- enabled-cell iteration order where public behavior relies on it;
- projection to selected dimensions;
- support-driven variable handling during projection;
- dimension duplication/removal behavior used by historical projection;
- variable shift/rename after projection.

Where the current Kernel tree implementation has evolved since `release-1.0`, the current `PavingInterface` semantics take precedence. Historical tests are an oracle for intended behavior, not a requirement to preserve obsolete implementation details.

## 16. Differential correctness strategy

Every operation that exists in both representations should be testable by applying the same logical workload to a tree paving and an ROBDD paving.

A useful pattern is:

```text
operation sequence
      |
      +-- GridTreePaving
      |
      +-- GridRobddPaving
              |
              v
       semantic comparison
```

Tests should include deterministic hand-built cases and randomized operation sequences.

After each operation, compare as applicable:

- enabled-cell set;
- `size()`;
- root-cell semantics;
- bounding box;
- subset in both directions;
- intersection/disjointness;
- branch results;
- iteration result;
- projection result;
- behavior after mince/recombine.

The tree backend should be treated as a reference implementation during ROBDD bring-up, but disagreements must still be classified. A failing comparison can be:

1. an ROBDD backend bug;
2. an existing tree backend bug;
3. an underspecified `PavingInterface` behavior;
4. a historical compatibility assumption that is no longer intended.

Tests must not be weakened merely to force representation agreement.

## 17. Benchmark strategy

Performance comparisons must separate high-level paving cost from low-level ROBDD cost.

For each representative workload record at least:

- wall-clock time;
- peak/resident memory where available;
- number of enabled logical cells;
- explicit tree node count;
- ROBDD non-terminal node count;
- ROBDD unique-table hits/misses;
- ROBDD computed-table hits/misses;
- operation count where useful;
- copy/clone cost;
- iteration cost;
- thread count and scaling for concurrent workloads.

Representative operations should include:

- construction by repeated cell insertion;
- bulk outer/inner approximation;
- union/adjoin;
- intersection/restrict;
- difference/remove;
- subset;
- intersects/disjoint;
- iteration;
- `size()`;
- branching/subpaving traversal;
- mince/recombine;
- projection;
- repeated snapshot/copy workloads.

Benchmark data should include both sparse and dense pavings, shallow and deep subdivisions, regular and irregular geometry, and cases with high subtree repetition where ROBDD reduction is expected to help.

The goal is not to prove that one representation universally dominates the other. The useful result may be a workload-dependent recommendation or an eventual storage-selection policy.

## 18. Concurrency comparison

The explicit tree implementation is currently mutable and pointer-based. The ROBDD representation naturally supports immutable canonical nodes and root replacement.

This difference should be measured rather than hidden.

ROBDD-oriented concurrent workloads should include:

- concurrent read-only queries against the same roots;
- concurrent independent operations sharing one manager;
- concurrent creation of equivalent nodes to stress canonicity;
- concurrent union/intersection/difference on independent paving roots;
- iterator/read workloads concurrent with creation of new roots.

The first ROBDD paving implementation should not introduce unsafe in-place mutation of BDD nodes merely to imitate tree semantics.

Higher-level `GridRobddPaving` remains mutable as a value/object by replacing its root after each successful operation.

## 19. API naming

The repository remains `ariadne-cps/bdd` because ROBDD and ZDD are binary decision-diagram families and the current generic layer is appropriately scoped.

Representation-specific classes should remain explicit:

```text
Robdd
RobddManager

GridRobddSubpaving
GridRobddPaving
```

Generic concepts such as `BddVariable` may remain shared in the BDD package.

The paving interfaces should not expose ROBDD terminology. Code that only needs geometric set semantics should continue to depend on `SubPavingInterface` or `PavingInterface`.

## 20. Suggested implementation phases

### Phase P0 - interface audit

- catalogue all mandatory casts from `SubPavingInterface`/`PavingInterface` to `GridTreeSubpaving`/`GridTreePaving`;
- catalogue all public/default aliases that hard-code `GridTreePaving`;
- document ambiguous `mince`, `recombine`, branch and iteration semantics;
- add cross-interface tests around the current tree implementation before changing behavior.

### Phase P1 - make the paving abstraction real

- retain tree/tree optimized operations;
- replace mandatory same-type casts with optional fast paths plus correct generic fallbacks;
- test operations through `PavingInterface&` and `SubPavingInterface&`;
- do not change the default paving representation yet.

### Phase P2 - complete the required ROBDD structural API

Implement and test in `ariadne-cps/bdd`:

- top-variable inspection;
- low/high cofactors;
- canonical decision-node construction or equivalent;
- support extraction;
- shift/rename/replacement;
- existential quantification;
- fused `and_exists`;
- efficient traversal with skipped variables;
- fixed-order manager policy;
- instrumentation needed by paving benchmarks.

### Phase P3 - read-only ROBDD subpaving

Implement `GridRobddSubpaving` first:

- root/grid metadata;
- branch views;
- iteration;
- `size()`;
- bounding box;
- cell containment;
- subset/intersection queries where practical.

Validate it against converted tree pavings.

### Phase P4 - mutable ROBDD paving

Implement `GridRobddPaving`:

- cell insertion/removal;
- union/adjoin;
- intersection/restrict;
- difference/remove;
- outer/inner approximation insertion;
- root-height alignment;
- mince/recombine semantics.

### Phase P5 - projection

Add historical projection semantics using:

- support;
- existential quantification;
- `and_exists`;
- variable shift/renaming.

Projection should receive dedicated differential tests because it exercises the most BDD-specific structural transformations.

### Phase P6 - performance and concurrency

- establish tree versus ROBDD benchmark suites;
- replace coarse BDD data structures only with measured motivation;
- add structured parallel recursion;
- run BDD and Kernel ThreadSanitizer stress tests;
- compare regular versus complemented edges for paving workloads.

### Phase P7 - choose defaults only from evidence

Only after correctness and benchmark data are stable should Kernel consider:

- changing `SetApproximationType`;
- adding an explicit paving-storage factory/policy;
- selecting a representation automatically by workload;
- keeping the tree backend as the permanent default while exposing ROBDD explicitly.

## 21. Open design questions

The following questions should remain explicit until implementation or benchmark evidence resolves them.

1. What is the ownership scope of `GridRobddContext` and its shared `RobddManager`?
2. Should all grids share a manager, or should manager/context identity be tied to a grid/dimension/order configuration?
3. What is the exact logical-depth metadata required to reproduce current mince/recombine and subpaving semantics?
4. What iterator order is part of the public contract, and must ROBDD iteration reproduce tree traversal order exactly?
5. What is the best generic fallback for heterogeneous `subset` and `intersects` without pathological cell expansion?
6. Should same-backend fast paths remain concrete overloads, use internal visitors, or eventually use double dispatch?
7. Which current `GridTreePaving`-specific aliases should remain compatibility defaults after the ROBDD backend exists?
8. When do complemented edges improve grid-paving workloads enough to justify enabling them in a paving manager?
9. What reclamation policy is needed for long-running paving computations with many obsolete roots?
10. Should projection remain a concrete `GridRobddPaving` operation initially, or be promoted into the common paving abstraction once both backends have a satisfactory implementation?

## 22. Immediate next steps

The immediate work should proceed in this order:

1. add tests exposing the current same-representation assumption in `PavingInterface`;
2. remove mandatory `GridTreeSubpaving` casts from generic operations while preserving tree/tree fast paths;
3. define the minimum read-only structural ROBDD API required by `GridRobddSubpaving`;
4. implement and test those ROBDD primitives in `ariadne-cps/bdd`;
5. introduce a minimal read-only `GridRobddSubpaving`;
6. build tree/ROBDD differential tests before adding the full mutating backend.

This order keeps geometry semantics and decision-diagram mechanics independently testable and avoids building a large ROBDD paving implementation on top of an interface that still assumes tree storage.
