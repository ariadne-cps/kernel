/***************************************************************************
 *            test_expression.cpp
 *
 *  Copyright  2009-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <sstream>
#include "utility/container.hpp"
#include "utility/stlio.hpp"
#include "numeric/numeric.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/templates.tpl.hpp"
#include "symbolic/expression.tpl.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/valuation.hpp"
#include "symbolic/space.hpp"
#include "function/formula.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/taylor_model.hpp"

#include "algebra/matrix.tpl.hpp"

#include "macro/testing.hpp"

using namespace Ariadne;

namespace Ariadne {
extern template class Expression<Real>;
}

typedef ElementaryAlgebra<Real> RealElementaryAlgebra;
typedef ElementaryAlgebra<EffectiveNumber> EffectiveElementaryAlgebra;
typedef SymbolicAlgebra<EffectiveNumber> EffectiveSymbolicAlgebra;

class TestExpression {
    RealConstant o;
    RealVariable x,y,z;
  public:
    TestExpression()
        : o("1.0",1.0_q), x("x"), y("y"), z("z") {
    }

    Void test_variables() {
        Identifier empty_identifier;
        ARIADNE_TEST_ASSERT(empty_identifier.empty());

        ARIADNE_TEST_CONSTRUCT(RealVariable,a,("a"));
        ARIADNE_TEST_ASSERT(a==RealVariable("a"));
        ARIADNE_TEST_ASSERT(a==RealVariable(a));
        ARIADNE_TEST_ASSERT(a!=RealVariable("b"));

        Variable<Boolean> boolean_variable("flag");
        Variable<Kleenean> kleenean_variable("predicate");
        ARIADNE_TEST_ASSERT(variable_type<Boolean>()==VariableType::BOOLEAN);
        ARIADNE_TEST_ASSERT(variable_type<Kleenean>()==VariableType::KLEENEAN);
        ARIADNE_TEST_ASSERT(variable_type<String>()==VariableType::STRING);
        ARIADNE_TEST_ASSERT(variable_type<Integer>()==VariableType::INTEGER);
        ARIADNE_TEST_ASSERT(variable_type<Real>()==VariableType::REAL);
        ARIADNE_TEST_ASSERT(variable_type<RealVector>()==VariableType::REAL_VECTOR);
        ARIADNE_TEST_ASSERT(variable_type<RealMatrix>()==VariableType::REAL_MATRIX);

        ARIADNE_TEST_EQUALS(class_name(VariableType::BOOLEAN),String("Boolean"));
        ARIADNE_TEST_EQUALS(class_name(VariableType::KLEENEAN),String("Kleenean"));
        ARIADNE_TEST_EQUALS(class_name(VariableType::ENUMERATED),String("Enumerated"));
        ARIADNE_TEST_EQUALS(class_name(VariableType::STRING),String("String"));
        ARIADNE_TEST_EQUALS(class_name(VariableType::INTEGER),String("Integer"));
        ARIADNE_TEST_EQUALS(class_name(VariableType::REAL),String("Real"));
        ARIADNE_TEST_EQUALS(class_name(VariableType::REAL_VECTOR),String("RealVector"));
        ARIADNE_TEST_EQUALS(class_name(VariableType::REAL_MATRIX),String("RealMatrix"));
        ARIADNE_TEST_FAIL((void)class_name(static_cast<VariableType>(127)));

        RealVariable real_same_name("same");
        Variable<Boolean> boolean_same_name("same");
        RealVariable real_a("a");
        RealVariable real_b("b");
        UntypedVariable const& untyped_real=real_same_name;
        UntypedVariable const& untyped_boolean=boolean_same_name;
        ARIADNE_TEST_ASSERT(untyped_real!=untyped_boolean);
        ARIADNE_TEST_ASSERT(static_cast<UntypedVariable const&>(real_a)<static_cast<UntypedVariable const&>(real_b));
        ARIADNE_TEST_ASSERT(untyped_boolean<untyped_real);

        ARIADNE_TEST_EQUALS(to_string(VariableCategory::SIMPLE),String("SIMPLE"));
        ARIADNE_TEST_EQUALS(to_string(VariableCategory::PRIMED),String("PRIMED"));
        ARIADNE_TEST_EQUALS(to_string(VariableCategory::DOTTED),String("DOTTED"));
        ARIADNE_TEST_FAIL((void)to_string(static_cast<VariableCategory>(127)));

        class InvalidExtendedVariable : public ExtendedUntypedVariable {
          public:
            InvalidExtendedVariable()
                : ExtendedUntypedVariable("invalid",VariableType::REAL,static_cast<VariableCategory>(127)) { }
        };
        ARIADNE_TEST_FAIL((void)to_string(InvalidExtendedVariable()));

        ARIADNE_TEST_EQUALS(to_string(let(x)),String("x"));
        ARIADNE_TEST_EQUALS(to_string(prime(x)),String("prime(x)"));
        ARIADNE_TEST_EQUALS(to_string(dot(x)),String("dot(x)"));

        List<RealVariable> variables={x,y,z};
        auto set_variables=set(variables);
        auto next_variables=next(variables);
        auto dotted_variables=dot(variables);
        ARIADNE_TEST_EQUALS(set_variables._lhs.size(),3u);
        ARIADNE_TEST_EQUALS(next_variables._lhs.size(),3u);
        ARIADNE_TEST_EQUALS(dotted_variables._lhs.size(),3u);

        auto set_initializer=set(InitializerList<RealVariable>{x,y});
        auto next_initializer=next(InitializerList<RealVariable>{x,y});
        ARIADNE_TEST_EQUALS(set_initializer._lhs.size(),2u);
        ARIADNE_TEST_EQUALS(next_initializer._lhs.size(),2u);

        ARIADNE_TEST_ASSERT(boolean_variable.name()==Identifier("flag"));
        ARIADNE_TEST_ASSERT(kleenean_variable.name()==Identifier("predicate"));
    }

    Void test_space() {
        RealSpace empty_space;
        ARIADNE_TEST_EQUALS(empty_space.size(),0u);

        List<RealVariable> variables_list={x,y};
        List<Identifier> names={Identifier("x"),Identifier("y")};
        RealSpace from_variables(variables_list);
        RealSpace from_names(names);
        RealSpace from_initializer({x,y});

        ARIADNE_TEST_ASSERT(from_variables==from_names);
        ARIADNE_TEST_ASSERT(not (from_variables!=from_names));
        ARIADNE_TEST_ASSERT(from_variables==from_initializer);
        RealSpace different({x,z});
        ARIADNE_TEST_ASSERT(from_variables!=different);
        ARIADNE_TEST_ASSERT(not (from_variables==different));

        auto name_indices=from_variables.indices_from_names();
        auto variable_indices=from_variables.indices();
        ARIADNE_TEST_EQUALS(name_indices[Identifier("x")],0u);
        ARIADNE_TEST_EQUALS(name_indices[Identifier("y")],1u);
        ARIADNE_TEST_EQUALS(variable_indices[x],0u);
        ARIADNE_TEST_EQUALS(variable_indices[y],1u);

        ARIADNE_TEST_ASSERT(from_variables.contains(x));
        ARIADNE_TEST_ASSERT(not from_variables.contains(z));
        Set<RealVariable> present;
        present.insert(x);
        present.insert(y);
        Set<RealVariable> missing=present;
        missing.insert(z);
        ARIADNE_TEST_ASSERT(from_variables.contains(present));
        ARIADNE_TEST_ASSERT(not from_variables.contains(missing));

        ARIADNE_TEST_EQUALS(from_variables[x],0u);
        ARIADNE_TEST_EQUALS(from_variables[Identifier("y")],1u);
        ARIADNE_TEST_EQUALS(from_variables.index(x),0u);
        ARIADNE_TEST_EQUALS(from_variables.index(Identifier("y")),1u);
        ARIADNE_TEST_FAIL((void)from_variables.index(z));
        ARIADNE_TEST_FAIL((void)from_variables.index(Identifier("z")));

        RealSpace inserted({x});
        inserted.insert(x);
        ARIADNE_TEST_EQUALS(inserted.size(),1u);
        inserted.insert(y);
        ARIADNE_TEST_EQUALS(inserted.size(),2u);

        RealSpace adjoined({x});
        adjoined.adjoin(RealSpace({x,y,z}));
        ARIADNE_TEST_EQUALS(adjoined.size(),3u);
        ARIADNE_TEST_ASSERT(adjoined.contains(z));

        RealSpace append_space({x});
        append_space.append(y);
        ARIADNE_TEST_EQUALS(append_space.size(),2u);
        ARIADNE_TEST_FAIL((void)append_space.append(y));
        ARIADNE_TEST_FAIL((void)RealSpace({x,x}));

        ARIADNE_TEST_EQUALS(dimension(from_variables),2u);
        ARIADNE_TEST_ASSERT(real_space(names)==from_variables);
        ARIADNE_TEST_ASSERT(variable_names(variables_list)==names);
        ARIADNE_TEST_ASSERT(variable_names(from_variables)==names);
    }

    Void test_constants() {
        StringConstant state(String("active"));
        ARIADNE_TEST_EQUALS(state.name(),Identifier("active"));
        ARIADNE_TEST_EQUALS(state.value(),String("active"));
        ARIADNE_TEST_EQUALS(state.val(),String("active"));
    }

    Void test_valuations() {
        Variable<String> mode("mode");
        Variable<String> state("state");
        Variable<Integer> count("count");
        Variable<Integer> level("level");

        Map<Identifier,String> string_values;
        string_values.insert(mode.name(),String("idle"));
        string_values.insert(state.name(),String("ready"));

        Map<Identifier,Integer> integer_values;
        integer_values.insert(count.name(),Integer(1));
        integer_values.insert(level.name(),Integer(2));

        DiscreteValuation empty_valuation;
        DiscreteValuation string_valuation(string_values);
        DiscreteValuation integer_valuation(integer_values);
        DiscreteValuation valuation(string_values,integer_values);

        ARIADNE_TEST_EQUALS(valuation[mode],String("idle"));
        ARIADNE_TEST_EQUALS(valuation[count],Integer(1));
        ARIADNE_TEST_EQUALS(valuation.string_values().size(),2u);
        ARIADNE_TEST_EQUALS(valuation.integer_values().size(),2u);

        StringValuation direct_string_valuation(string_values);
        IntegerValuation direct_integer_valuation(integer_values);

        std::ostringstream string_stream;
        string_stream << direct_string_valuation;
        ARIADNE_TEST_ASSERT(not string_stream.str().empty());

        std::ostringstream integer_stream;
        integer_stream << direct_integer_valuation;
        ARIADNE_TEST_ASSERT(not integer_stream.str().empty());

        std::ostringstream discrete_stream;
        discrete_stream << valuation;
        ARIADNE_TEST_ASSERT(not discrete_stream.str().empty());

        ARIADNE_TEST_ASSERT(direct_string_valuation==StringValuation(string_values));

        Map<Identifier,String> different_string_values=string_values;
        different_string_values[mode.name()]=String("busy");
        ARIADNE_TEST_ASSERT(not (direct_string_valuation==StringValuation(different_string_values)));

        Map<Identifier,String> lower_key_string_values;
        lower_key_string_values.insert(Identifier("a"),String("same"));
        Map<Identifier,String> upper_key_string_values;
        upper_key_string_values.insert(Identifier("b"),String("same"));
        ARIADNE_TEST_FAIL((void)(StringValuation(lower_key_string_values)==StringValuation(upper_key_string_values)));
        ARIADNE_TEST_FAIL((void)(StringValuation(upper_key_string_values)==StringValuation(lower_key_string_values)));

        DiscreteValuation equal_valuation(string_values,integer_values);
        ARIADNE_TEST_ASSERT(valuation==equal_valuation);

        DiscreteValuation different_string_valuation(different_string_values,integer_values);
        ARIADNE_TEST_ASSERT(not (valuation==different_string_valuation));

        Map<Identifier,Integer> different_integer_values=integer_values;
        different_integer_values[count.name()]=Integer(3);
        DiscreteValuation different_integer_valuation(string_values,different_integer_values);
        ARIADNE_TEST_ASSERT(not (valuation==different_integer_valuation));

        ARIADNE_TEST_ASSERT(empty_valuation.string_values().empty());
        ARIADNE_TEST_ASSERT(empty_valuation.integer_values().empty());
        ARIADNE_TEST_ASSERT(string_valuation.integer_values().empty());
        ARIADNE_TEST_ASSERT(integer_valuation.string_values().empty());
    }

    Void test_templates() {
        Get get;
        ARIADNE_TEST_EQUALS(get.kind(),OperatorKind::COORDINATE);

        Vector<Int> values={3,5};
        ARIADNE_TEST_EQUALS(get(values,1u),Int(5));

        Vec vec;
        ARIADNE_TEST_EQUALS(vec.code(),OperatorCode::VEC);
        ARIADNE_TEST_EQUALS(vec.kind(),OperatorKind::NULLARY);
        ARIADNE_TEST_ASSERT(not to_string(vec).empty());

        Symbolic<Add,Int,Int> binary(Add(),Int(1),Int(2));
        std::ostringstream binary_stream;
        binary_stream << binary;
        ARIADNE_TEST_ASSERT(not binary_stream.str().empty());

        Symbolic<Get,Vector<Int>,SizeType> getter(Get(),values,1u);
        ARIADNE_TEST_EQUALS(getter.op().code(),OperatorCode::GET);
        ARIADNE_TEST_EQUALS(getter.arg().size(),2u);
        ARIADNE_TEST_EQUALS(getter.arg1().size(),2u);
        ARIADNE_TEST_EQUALS(getter.arg2(),1u);
        ARIADNE_TEST_EQUALS(getter.vec().size(),2u);
        ARIADNE_TEST_EQUALS(getter.ind(),1u);

        RealVariables template_variables("template",2u);
        RealVectorExpression vector_expression=template_variables;
        RealExpression direct_component=get(vector_expression,1u);

        using GetMethod=RealExpression (Get::*)(RealVectorExpression const&,SizeType);
        GetMethod get_method=&Get::operator()<RealVectorExpression,SizeType>;
        RealExpression functor_component=(get.*get_method)(vector_expression,SizeType(1));
        ARIADNE_TEST_ASSERT(identical(functor_component,direct_component));

        using GetterOperator=OperatorVariant<Get>;
        using GetterNode=Symbolic<GetterOperator,RealVectorExpression,SizeType>;
        GetterNode concrete_getter(GetterOperator(Get()),vector_expression,SizeType(1));

        using GetterOpMethod=GetterOperator (GetterNode::*)() const;
        using GetterVecMethod=RealVectorExpression (GetterNode::*)() const;
        using GetterIndMethod=SizeType (GetterNode::*)() const;
        GetterOpMethod getter_op_method=&GetterNode::op;
        GetterVecMethod getter_vec_method=&GetterNode::vec;
        GetterIndMethod getter_ind_method=&GetterNode::ind;

        ARIADNE_TEST_EQUALS((concrete_getter.*getter_op_method)().code(),OperatorCode::GET);
        ARIADNE_TEST_EQUALS((concrete_getter.*getter_vec_method)().size(),2u);
        ARIADNE_TEST_EQUALS((concrete_getter.*getter_ind_method)(),SizeType(1));

    }

    Void test_scalar_expression_templates() {
        CatOp cat;
        ARIADNE_TEST_EQUALS(cat.code(),OperatorCode::ADD);
        ARIADNE_TEST_EQUALS(cat.kind(),OperatorKind::BINARY);
        ARIADNE_TEST_EQUALS(cat(String("a"),String("b")),String("ab"));
        ARIADNE_TEST_EQUALS(to_string(cat),String("cat"));
        ARIADNE_TEST_EQUALS(pos(String("identity")),String("identity"));

        BooleanVariable boolean_variable("boolean");
        BooleanVariable boolean_other("boolean_other");
        BooleanExpression boolean_default;
        BooleanExpression boolean_copy(boolean_default);
        BooleanExpression boolean_constant{Constant<Boolean>(Boolean())};
        ARIADNE_TEST_ASSERT(identical(boolean_default,boolean_copy));
        ARIADNE_TEST_EXECUTE(boolean_constant.val());
        BooleanExpression boolean_from_variable(boolean_variable);
        ARIADNE_TEST_EXECUTE(boolean_default.val());
        ARIADNE_TEST_EQUALS(boolean_from_variable.var(),Identifier("boolean"));
        ARIADNE_TEST_EXECUTE(boolean_from_variable.op());
        ARIADNE_TEST_EXECUTE(boolean_from_variable.code());
        ARIADNE_TEST_EXECUTE(boolean_from_variable.kind());
        ARIADNE_TEST_EQUALS(boolean_from_variable.arguments().size(),1u);
        ARIADNE_TEST_EQUALS(arguments(boolean_from_variable).size(),1u);
        ARIADNE_TEST_EXECUTE(BooleanExpression::constant(Boolean()));
        ARIADNE_TEST_EXECUTE(BooleanExpression::constant(Constant<Boolean>(Boolean())));
        ARIADNE_TEST_EXECUTE(BooleanExpression::variable(Identifier("boolean")));
        BooleanExpression boolean_not=!boolean_from_variable;
        BooleanExpression boolean_binary=boolean_from_variable&&BooleanExpression(boolean_other);
        ARIADNE_TEST_EQUALS(boolean_not.arg().var(),Identifier("boolean"));
        ARIADNE_TEST_EQUALS(boolean_binary.arg1().var(),Identifier("boolean"));
        ARIADNE_TEST_EQUALS(boolean_binary.arg2().var(),Identifier("boolean_other"));
        ARIADNE_TEST_ASSERT(not to_string(boolean_binary).empty());

        KleeneanVariable kleenean_variable("kleenean");
        KleeneanVariable kleenean_other("kleenean_other");
        KleeneanExpression kleenean_default;
        KleeneanExpression kleenean_copy(kleenean_default);
        KleeneanExpression kleenean_constant{Constant<Kleenean>(Kleenean())};
        ARIADNE_TEST_ASSERT(identical(kleenean_default,kleenean_copy));
        ARIADNE_TEST_EXECUTE(kleenean_constant.val());
        KleeneanExpression kleenean_from_variable(kleenean_variable);
        ARIADNE_TEST_EXECUTE(kleenean_default.val());
        ARIADNE_TEST_EQUALS(kleenean_from_variable.var(),Identifier("kleenean"));
        ARIADNE_TEST_EXECUTE(kleenean_from_variable.op());
        ARIADNE_TEST_EXECUTE(kleenean_from_variable.code());
        ARIADNE_TEST_EXECUTE(kleenean_from_variable.kind());
        ARIADNE_TEST_EQUALS(kleenean_from_variable.arguments().size(),1u);
        ARIADNE_TEST_EXECUTE(KleeneanExpression::constant(Kleenean()));
        ARIADNE_TEST_EXECUTE(KleeneanExpression::constant(Constant<Kleenean>(Kleenean())));
        ARIADNE_TEST_EXECUTE(KleeneanExpression::variable(Identifier("kleenean")));
        KleeneanExpression kleenean_not=!kleenean_from_variable;
        KleeneanExpression kleenean_binary=kleenean_from_variable&&KleeneanExpression(kleenean_other);
        ARIADNE_TEST_EQUALS(kleenean_not.arg().var(),Identifier("kleenean"));
        ARIADNE_TEST_EQUALS(kleenean_binary.arg1().var(),Identifier("kleenean"));
        ARIADNE_TEST_EQUALS(kleenean_binary.arg2().var(),Identifier("kleenean_other"));
        ARIADNE_TEST_ASSERT(same(Kleenean(),Kleenean()));
        ARIADNE_TEST_ASSERT(not to_string(kleenean_binary).empty());

        StringVariable string_variable("string");
        StringExpression string_default;
        StringExpression string_value(String("value"));
        StringExpression string_copy(string_value);
        ARIADNE_TEST_EQUALS(string_default.val(),String());
        ARIADNE_TEST_ASSERT(identical(string_value,string_copy));
        StringExpression string_constant(StringConstant(String("named")));
        StringExpression string_from_variable(string_variable);
        ARIADNE_TEST_EQUALS(string_value.val(),String("value"));
        ARIADNE_TEST_EQUALS(string_from_variable.var(),Identifier("string"));
        ARIADNE_TEST_EXECUTE(string_from_variable.op());
        ARIADNE_TEST_EXECUTE(string_from_variable.code());
        ARIADNE_TEST_EXECUTE(string_from_variable.kind());
        ARIADNE_TEST_EQUALS(string_from_variable.arguments().size(),1u);
        ARIADNE_TEST_EXECUTE(StringExpression::constant(String("value")));
        ARIADNE_TEST_EXECUTE(StringExpression::constant(StringConstant(String("named"))));
        ARIADNE_TEST_EXECUTE(StringExpression::variable(Identifier("string")));
        ARIADNE_TEST_ASSERT(not to_string(string_from_variable).empty());

        IntegerVariable integer_variable("integer");
        IntegerVariable integer_other("integer_other");
        IntegerExpression integer_default;
        IntegerExpression integer_copy(integer_default);
        ARIADNE_TEST_ASSERT(identical(integer_default,integer_copy));
        IntegerExpression integer_constant(IntegerConstant(Integer(7)));
        IntegerExpression integer_from_variable(integer_variable);
        ARIADNE_TEST_EQUALS(integer_constant.val(),Integer(7));
        ARIADNE_TEST_EQUALS(integer_from_variable.var(),Identifier("integer"));
        ARIADNE_TEST_EXECUTE(integer_from_variable.op());
        ARIADNE_TEST_EXECUTE(integer_from_variable.code());
        ARIADNE_TEST_EXECUTE(integer_from_variable.kind());
        ARIADNE_TEST_EQUALS(integer_from_variable.arguments().size(),1u);
        ARIADNE_TEST_EXECUTE(IntegerExpression::constant(Integer(7)));
        ARIADNE_TEST_EXECUTE(IntegerExpression::constant(IntegerConstant(Integer(7))));
        ARIADNE_TEST_EXECUTE(IntegerExpression::variable(Identifier("integer")));
        IntegerExpression integer_negative=-integer_from_variable;
        IntegerExpression integer_binary=integer_from_variable+IntegerExpression(integer_other);
        ARIADNE_TEST_EQUALS(integer_negative.arg().var(),Identifier("integer"));
        ARIADNE_TEST_EQUALS(integer_binary.arg1().var(),Identifier("integer"));
        ARIADNE_TEST_EQUALS(integer_binary.arg2().var(),Identifier("integer_other"));
        ARIADNE_TEST_ASSERT(not to_string(integer_binary).empty());

        RealExpression real_constant(Real(2));
        ARIADNE_TEST_EQUAL(real_constant.val(),Real(2));
        RealExpression real_power=pow(x,3);
        ARIADNE_TEST_EQUALS(real_power.num(),3);

        StringValuation string_valuation;
        string_valuation.set(string_variable,String("mapped"));
        ARIADNE_TEST_EQUALS(evaluate(string_from_variable,string_valuation),String("mapped"));
        ARIADNE_TEST_EXECUTE(evaluate(string_variable==String("mapped"),string_valuation));

        IntegerValuation integer_valuation;
        integer_valuation.set(integer_variable,Integer(9));
        ARIADNE_TEST_EQUALS(evaluate(integer_from_variable,integer_valuation),Integer(9));
        ARIADNE_TEST_EXECUTE(evaluate(integer_from_variable<IntegerExpression(integer_other),integer_valuation));

        RealValuation real_valuation;
        real_valuation.set(x,Real(2));
        ARIADNE_TEST_EXECUTE(evaluate(x>=0,real_valuation));

        Map<Identifier,String> string_map;
        string_map[string_variable.name()]=String("mapped");
        ARIADNE_TEST_EQUALS(evaluate(string_from_variable,string_map),String("mapped"));
        ARIADNE_TEST_EQUALS(evaluate(StringExpression(String("fixed")),string_map),String("fixed"));
        ARIADNE_TEST_EQUALS(evaluate(IntegerExpression(Integer(5)),string_map),Integer(5));

        Map<Identifier,Integer> integer_map;
        integer_map[integer_variable.name()]=Integer(9);
        ARIADNE_TEST_EQUALS(evaluate(integer_from_variable,integer_map),Integer(9));
        ARIADNE_TEST_EQUALS(evaluate(IntegerExpression(Integer(5)),integer_map),Integer(5));
        ARIADNE_TEST_EQUALS(evaluate(StringExpression(String("fixed")),integer_map),String("fixed"));
    }

    Void test_expression() {
        auto writer=RealExpression::default_writer();
        RealExpression::set_default_writer(writer);

        // Regression test for constructing Expression from 0 without being an ambiguous nullptr;
        ARIADNE_TEST_CONSTRUCT(IntegerExpression,ze,(0));
        ARIADNE_TEST_CONSTRUCT(RealExpression,re,(0));
        RealExpression(0);
        RealExpression(nullptr);
    }

    Void test_write() {
        ARIADNE_TEST_EQUALS(to_string(+x),"+x");
        ARIADNE_TEST_EQUALS(to_string(-x),"-x");
        ARIADNE_TEST_EQUALS(to_string(x+y),"x+y");
        ARIADNE_TEST_EQUALS(to_string(x-y),"x-y");
        ARIADNE_TEST_EQUALS(to_string(x*y),"x*y");
        ARIADNE_TEST_EQUALS(to_string(x/y),"x/y");
        ARIADNE_TEST_EQUALS(to_string(pos(x)),"+x");
        ARIADNE_TEST_EQUALS(to_string(neg(x)),"-x");
        ARIADNE_TEST_EQUALS(to_string(add(x,y)),"x+y");
        ARIADNE_TEST_EQUALS(to_string(sub(x,y)),"x-y");
        ARIADNE_TEST_EQUALS(to_string(mul(x,y)),"x*y");
        ARIADNE_TEST_EQUALS(to_string(div(x,y)),"x/y");
        ARIADNE_TEST_EQUALS(to_string(neg(sub(x,y))),"-(x-y)");
        ARIADNE_TEST_EQUALS(to_string(sub(neg(x),y)),"-x-y");
        ARIADNE_TEST_EQUALS(to_string(add(x,sub(y,z))),"x+y-z");
        ARIADNE_TEST_EQUALS(to_string(sub(x,add(y,z))),"x-(y+z)");
        ARIADNE_TEST_EQUALS(to_string(sub(x,sub(y,z))),"x-(y-z)");
        ARIADNE_TEST_EQUALS(to_string(add(sub(x,y),z)),"x-y+z");
        ARIADNE_TEST_EQUALS(to_string(sub(add(x,y),z)),"x+y-z");
        ARIADNE_TEST_EQUALS(to_string(sub(sub(x,y),z)),"x-y-z");
        ARIADNE_TEST_EQUALS(to_string(add(mul(x,y),z)),"x*y+z");
        ARIADNE_TEST_EQUALS(to_string(mul(add(x,y),z)),"(x+y)*z");
        ARIADNE_TEST_EQUALS(to_string(div(x,div(y,z))),"x/(y/z)");
        ARIADNE_TEST_EQUALS(to_string(div(div(x,y),z)),"x/y/z");
        ARIADNE_TEST_EQUALS(to_string(sub(x,div(y,z))),"x-y/z");
        ARIADNE_TEST_EQUALS(to_string(sub(div(x,y),z)),"x/y-z");
        ARIADNE_TEST_EQUALS(to_string(max(x,sub(y,z))),"max(x,y-z)");
        ARIADNE_TEST_EQUALS(to_string(sub(x,max(y,z))),"x-max(y,z)");
        ARIADNE_TEST_EQUALS(to_string(tanh(x)),"tanh(x)");
    }

    Void test_assignment() {
        Real zero(0), one(1);
        RealExpression e(x*y+o);

        ARIADNE_TEST_ASSERT((not Assignable<Variable<Real>,Expression<Real>>));
        ARIADNE_TEST_ASSERT((not Constructible<Assignment<Variable<Real>,Expression<Real>>,Assignment<Variable<Real>,Real>>));

        typedef Assignment<Variable<Real>,Real> ConstantRealAssignment;
        ARIADNE_TEST_CONSTRUCT(ConstantRealAssignment,ac,(x=one));
        ARIADNE_TEST_CONSTRUCT(List<ConstantRealAssignment>,lac,({x=zero,y=one}));
        ARIADNE_TEST_CONSTRUCT(Valuation<Real>,va,(lac));
        //ARIADNE_TEST_CONSTRUCT(Valuation<Real>,va,({x=zero,y=one})); Fails due to ambiguous overload

        ARIADNE_TEST_CONSTRUCT(RealAssignment,a,(let(x)=one));
        ARIADNE_TEST_CONSTRUCT(PrimedRealAssignment,pa,(prime(x)=one));
        ARIADNE_TEST_CONSTRUCT(DottedRealAssignment,da,(dot(x)=one));

        ARIADNE_TEST_CONSTRUCT(List<RealAssignment>,la,(let({x,y,z})={zero,x,e}));
        ARIADNE_TEST_CONSTRUCT(List<PrimedRealAssignment>,lpa,(prime({x,y,z})={zero,x,e}));
        ARIADNE_TEST_CONSTRUCT(List<DottedRealAssignment>,lda,(dot({x,y,z})={zero,x,e}));

        FloatDPApproximation approximate_value(1.0_x,dp);
        Assignment<RealVariable,FloatDPApproximation> approximate_assignment(x,approximate_value);
        ARIADNE_TEST_ASSERT(&approximate_assignment.right_hand_side()==&approximate_assignment.rhs);
    }

    Void test_evaluate() {
        ARIADNE_TEST_CONSTRUCT(RealExpression,g,(x+3*y*z*z));

        Map<RealVariable,Real> v;
        v[x]=Real(2.0_x); v[y]=Real(3.0_x); v[z]=Real(5.0_x);

        ARIADNE_TEST_PRINT(v);
        //ARIADNE_TEST_EQUAL(evaluate(g,v),Real(227));
    }

    Void test_parameters() {
        RealExpression expr = x;//+u;

        Map<Identifier,Real> valuation;
        Real value = Real(ExactDouble(-0.0626));
        valuation[x.name()] = value;

        ARIADNE_TEST_EQUALS(expr.kind(),OperatorKind::VARIABLE);
        ARIADNE_TEST_EQUALS(expr.var(),"x");
        ARIADNE_TEST_EQUALS(valuation[x.name()],value);

        Real result1 = evaluate(expr,valuation);

        ARIADNE_TEST_EQUALS(result1,value);
    }

    Void test_identical() {
        ARIADNE_TEST_ASSERT(identical(RealExpression(x),RealExpression(x)));
        ARIADNE_TEST_ASSERT(identical(RealExpression::constant(0),RealExpression::constant(0)));
        ARIADNE_TEST_ASSERT(identical(RealExpression::constant(2),RealExpression::constant(2)));
        ARIADNE_TEST_ASSERT(identical(sin(x),sin(x)));
        ARIADNE_TEST_ASSERT(identical(pow(x,2),pow(x,2)));
        ARIADNE_TEST_ASSERT(not identical(x*y,y*x));
        ARIADNE_TEST_ASSERT(not identical(x+y,y+x));

        ARIADNE_TEST_CONSTRUCT(RealExpression,e1,(sin(pow(x,2)+y)*(y/x+1)))
        ARIADNE_TEST_CONSTRUCT(RealExpression,e2,(sin(pow(x,2)+y)*(y/x+1)));
        ARIADNE_TEST_ASSERT(identical(e1,e1));
        ARIADNE_TEST_ASSERT(identical(e1,e2));
    }

    Void test_derivative() {
        RealExpression expr = 2*x+y;
        ARIADNE_TEST_PRINT(derivative(expr,x));
        ARIADNE_TEST_PRINT(simplify(derivative(expr,x)));
        ARIADNE_TEST_ASSERT(identical(simplify(derivative(expr,x)),RealExpression::constant(2)));
        RealExpression tanh_expr=tanh(x);
        ARIADNE_TEST_ASSERT(identical(
            simplify(derivative(tanh_expr,x)),
            1-sqr(tanh(x))));
        RealExpression expr2 = pow(x,3);
        ARIADNE_TEST_PRINT(derivative(expr2,x));
        ARIADNE_TEST_PRINT(simplify(derivative(expr2,x)));
        ARIADNE_TEST_ASSERT(identical(simplify(derivative(expr2,x)),3*sqr(x)));
    }

    Void test_simplify() {
        RealExpression zero(0);
        RealExpression one(1);
        RealVariable u("u");
        RealExpression expr = -u*x*y+2*x;
        RealExpression ex=x;
        ARIADNE_TEST_BINARY_PREDICATE(identical,simplify(derivative(-u*x*y+2*x,x)),-u*y+2);
        ARIADNE_TEST_BINARY_PREDICATE(identical,simplify(pos(neg(ex))),neg(ex));
        ARIADNE_TEST_BINARY_PREDICATE(identical,simplify(neg(neg(ex))),ex);
        ARIADNE_TEST_BINARY_PREDICATE(identical,simplify(rec(rec(ex))),ex);
        ARIADNE_TEST_BINARY_PREDICATE(identical,simplify(log(exp(ex))),ex);
        ARIADNE_TEST_BINARY_PREDICATE(identical,simplify(exp(log(ex))),ex);
        ARIADNE_TEST_BINARY_PREDICATE(identical,simplify(sqr(sqrt(ex))),ex);
        // Regression tests
        ARIADNE_TEST_BINARY_PREDICATE(identical,simplify(sub(neg(zero),zero)),zero);
        ARIADNE_TEST_BINARY_PREDICATE(identical,simplify(sub(x,x)),zero);
        ARIADNE_TEST_BINARY_PREDICATE(identical,simplify(div(x,x)),one);
    }

    Void test_ordering() {
        ARIADNE_TEST_ASSERT(before(RealExpression(x),RealExpression(y)));
        ARIADNE_TEST_ASSERT(not before(RealExpression(x),RealExpression(x)));
        ARIADNE_TEST_ASSERT(before(RealExpression(1),RealExpression(2)));
        ARIADNE_TEST_ASSERT(not before(RealExpression(5),RealExpression(2)));
        ARIADNE_TEST_ASSERT(before(pow(x,2),pow(x,3)));
        ARIADNE_TEST_ASSERT(before(pow(y,2),pow(x,3)));
        ARIADNE_TEST_ASSERT(before(pow(x,2),pow(x,3)));
        ARIADNE_TEST_ASSERT(not before(sin(y),sin(x)));
        ARIADNE_TEST_ASSERT(before(sqrt(x),sin(y)));
        ARIADNE_TEST_ASSERT(before(x+y,y+x));
        ARIADNE_TEST_ASSERT(not before(x+y,x+x));
        ARIADNE_TEST_ASSERT(not before(log(x),rec(x)));
    }

    Void test_count_nodes() {
        ARIADNE_TEST_EQUALS(count_nodes(RealExpression(x)),1);
        ARIADNE_TEST_EQUALS(count_nodes(RealExpression(1)),1);
        ARIADNE_TEST_EQUALS(count_nodes(2*x),3);
        ARIADNE_TEST_EQUALS(count_nodes(x+y),3);
        ARIADNE_TEST_EQUALS(count_nodes(x+sqr(y)),4);
        ARIADNE_TEST_EQUALS(count_nodes(sin(y)),2);
        ARIADNE_TEST_EQUALS(count_nodes(pow(y,2)),2);
        ARIADNE_TEST_EQUALS(count_nodes(pow(x+cos(y),2)+cos(y)),8);
        ARIADNE_TEST_EQUALS(count_nodes(x+cos(x)+pow(cos(x),2)),8);
    }

    Void test_count_distinct_nodes() {
        ARIADNE_TEST_EQUALS(count_distinct_nodes(RealExpression(x)),1);
        ARIADNE_TEST_EQUALS(count_distinct_nodes(x*exp(x)),3);
        ARIADNE_TEST_EQUALS(count_distinct_nodes(x*x),2);
        RealExpression one = 1;
        ARIADNE_TEST_EQUALS(count_distinct_nodes(one+cos(one)),3);
        ARIADNE_TEST_EQUALS(count_distinct_nodes(x*y+sqr(x*y)),5);
        ARIADNE_TEST_EQUALS(count_distinct_nodes(pow(x,2)*pow(x,2)),3);
        ARIADNE_TEST_EQUALS(count_distinct_nodes(pow(y,2)+y*y),4);
        ARIADNE_TEST_EQUALS(count_distinct_nodes(pow(x+cos(y),2)+cos(y)),6);
        ARIADNE_TEST_EQUALS(count_distinct_nodes(x+cos(x)+pow(cos(x),2)),5);
    }

    Void test_count_distinct_node_pointers() {
        ARIADNE_TEST_EQUALS(count_distinct_node_pointers(RealExpression(x)),1);
        ARIADNE_TEST_EQUALS(count_distinct_node_pointers(x*exp(x)),4);
        ARIADNE_TEST_EQUALS(count_distinct_node_pointers(x*x),3);
        RealExpression one = 1;
        ARIADNE_TEST_EQUALS(count_distinct_node_pointers(one+cos(one)),3);
        ARIADNE_TEST_EQUALS(count_distinct_node_pointers(x*y+sqr(x*y)),8);
        ARIADNE_TEST_EQUALS(count_distinct_node_pointers(pow(x,2)*pow(x,2)),5);
        ARIADNE_TEST_EQUALS(count_distinct_node_pointers(pow(y,2)+y*y),6);
        ARIADNE_TEST_EQUALS(count_distinct_node_pointers(pow(x+cos(y),2)+cos(y)),8);
        ARIADNE_TEST_EQUALS(count_distinct_node_pointers(x+cos(x)+pow(cos(x),2)),8);
    }

    Void test_eliminate_common_subexpressions() {
        RealExpression expr1 = x;
        ARIADNE_TEST_PRINT(expr1);
        eliminate_common_subexpressions(expr1);
        ARIADNE_TEST_EQUAL(count_distinct_node_pointers(expr1),1);
        RealExpression expr2 = x*exp(x);
        ARIADNE_TEST_PRINT(expr2);
        eliminate_common_subexpressions(expr2);
        ARIADNE_TEST_EQUAL(count_distinct_node_pointers(expr2),3);
        RealExpression expr3 = x*x;
        ARIADNE_TEST_PRINT(expr3);
        eliminate_common_subexpressions(expr3);
        eliminate_common_subexpressions(expr3);
        ARIADNE_TEST_PRINT(expr3);
        ARIADNE_TEST_EQUAL(count_distinct_node_pointers(expr3),2);
        RealExpression one = 1;
        RealExpression expr4 = one+cos(one);
        ARIADNE_TEST_PRINT(expr4);
        eliminate_common_subexpressions(expr4);
        ARIADNE_TEST_EQUAL(count_distinct_node_pointers(expr4),3);
        RealExpression expr5 = x*y+sqr(x*y);
        ARIADNE_TEST_PRINT(expr5);
        eliminate_common_subexpressions(expr5);
        ARIADNE_TEST_EQUAL(count_distinct_node_pointers(expr5),5);
        RealExpression expr6 = pow(x,2)*pow(x,2);
        ARIADNE_TEST_PRINT(expr6);
        eliminate_common_subexpressions(expr6);
        ARIADNE_TEST_EQUAL(count_distinct_node_pointers(expr6),3);
        RealExpression expr7 = pow(y,2)+y*y;
        ARIADNE_TEST_PRINT(expr7);
        eliminate_common_subexpressions(expr7);
        ARIADNE_TEST_EQUAL(count_distinct_node_pointers(expr7),4);
        RealExpression expr8 = pow(x+cos(y),2)+cos(y);
        ARIADNE_TEST_PRINT(expr8);
        eliminate_common_subexpressions(expr8);
        ARIADNE_TEST_EQUAL(count_distinct_node_pointers(expr8),6);
        RealExpression expr9 = x+cos(x)+pow(cos(x),2);
        ARIADNE_TEST_PRINT(expr9);
        eliminate_common_subexpressions(expr9);
        ARIADNE_TEST_EQUAL(count_distinct_node_pointers(expr9),5);

        // Cross-output CSE must survive Expression -> Formula conversion.
        Vector<RealExpression> expressions({x+y,x+y});
        ARIADNE_TEST_ASSERT(
            expressions[0].node_raw_ptr()!=expressions[1].node_raw_ptr());
        eliminate_common_subexpressions(expressions);
        ARIADNE_TEST_ASSERT(
            expressions[0].node_raw_ptr()==expressions[1].node_raw_ptr());
        auto first_root=expressions[0].node_raw_ptr();
        eliminate_common_subexpressions(expressions);
        ARIADNE_TEST_ASSERT(expressions[0].node_raw_ptr()==first_root);
        ARIADNE_TEST_ASSERT(
            expressions[0].node_raw_ptr()==expressions[1].node_raw_ptr());
        Vector<Formula<EffectiveNumber>> formulae=
            make_formula(expressions,RealSpace({x,y,z}));
        ARIADNE_TEST_ASSERT(
            formulae[0].node_ptr()==formulae[1].node_ptr());
    }

    Void test_substitute() {

        RealVariable u1("u1"), u2("u2");
        RealExpression expr = -u1*x*y+2*pow(x+u2,2);

        List<Assignment<RealVariable,RealExpression>> subs={{u1,u1+1},{u2,u1*x}};

        RealExpression substitution = substitute(expr,subs);

        ARIADNE_TEST_ASSERT(identical(substitution,-(u1+1)*x*y+2*pow(x+u1*x,2)));
    }

    Void test_is_constant_in() {
        Real c(3);
        ARIADNE_TEST_ASSERT(is_constant_in(3*y,{x}));
        ARIADNE_TEST_ASSERT(is_constant_in(pow(x,2),{y}));
        ARIADNE_TEST_ASSERT(not is_constant_in(pow(x,2),{x}));
        ARIADNE_TEST_ASSERT(not is_constant_in(3*y,{y}));
        ARIADNE_TEST_ASSERT(not is_constant_in(0*y,{y}));
        ARIADNE_TEST_ASSERT(not is_constant_in((sin(2*c)-2*sin(c)*cos(c))*y,{y}));
        ARIADNE_TEST_ASSERT(not is_constant_in((sin(2*x)-2*sin(x)*cos(x))*y,{y}));
        ARIADNE_TEST_ASSERT(is_constant_in(simplify(0*y),{y}));
    }

    Void test_is_additive_in() {
        RealVariable u1("u1"), u2("u2");
        ARIADNE_TEST_ASSERT(is_additive_in(u1,u1));
        ARIADNE_TEST_ASSERT(is_additive_in(x,u1));
        ARIADNE_TEST_ASSERT(is_additive_in(x+u1,u1));
        ARIADNE_TEST_ASSERT(is_additive_in(x+u1,u2));
        ARIADNE_TEST_ASSERT(is_additive_in(Vector<RealExpression>({x+u1}),{u1}));
        ARIADNE_TEST_ASSERT(is_additive_in(Vector<RealExpression>({x+u1}),{u2}));
        ARIADNE_TEST_ASSERT(is_additive_in(Vector<RealExpression>({x+u1,y+u2}),{u1,u2}));
        ARIADNE_TEST_ASSERT(is_additive_in(Vector<RealExpression>({x+u2,y+u1}),{u1,u2}));
        ARIADNE_TEST_ASSERT(is_additive_in(Vector<RealExpression>({x+u1,y}),{u1}));
        ARIADNE_TEST_ASSERT(is_additive_in(Vector<RealExpression>({x,y+u1}),{u1}));
        ARIADNE_TEST_ASSERT(not is_additive_in(Vector<RealExpression>({x,y+2*u1}),{u1}));
        ARIADNE_TEST_ASSERT(not is_additive_in(Vector<RealExpression>({x+u1,y+2*u2}),{u1,u2}));
        ARIADNE_TEST_ASSERT(not is_additive_in(Vector<RealExpression>({x+u1,y+u1}),{u1}));
        ARIADNE_TEST_ASSERT(not is_additive_in(Vector<RealExpression>({x*u1,y+u2}),{u1,u2}));
        ARIADNE_TEST_ASSERT(not is_additive_in(Vector<RealExpression>({x+u1,y+sqr(u2)}),{u1,u2}));
    }

    Void test_is_affine_in() {
        ARIADNE_TEST_ASSERT(is_affine_in(sqr(x),{y}));
        ARIADNE_TEST_ASSERT(is_affine_in(neg(x),{x}));
        ARIADNE_TEST_ASSERT(is_affine_in(pow(x,3),{y}));
        ARIADNE_TEST_ASSERT(is_affine_in(pow(x,3)+y,{y}));
        ARIADNE_TEST_ASSERT(is_affine_in(2+3*x-5*y-x,{x,y}));
        ARIADNE_TEST_ASSERT(is_affine_in(3*y,{x,y}));
        ARIADNE_TEST_ASSERT(is_affine_in(x*y,{x}));
        ARIADNE_TEST_ASSERT(is_affine_in(3*x/y,{x}));
        ARIADNE_TEST_ASSERT(not is_affine_in(atan(x),{x}))
        ARIADNE_TEST_ASSERT(not is_affine_in(acos(x),{x}))
        ARIADNE_TEST_ASSERT(not is_affine_in(asin(x),{x}))
        ARIADNE_TEST_ASSERT(not is_affine_in(cos(x),{x}))
        ARIADNE_TEST_ASSERT(not is_affine_in(sin(x),{x}))
        ARIADNE_TEST_ASSERT(not is_affine_in(tan(x),{x}))
        ARIADNE_TEST_ASSERT(not is_affine_in(abs(x),{x}))
        ARIADNE_TEST_ASSERT(not is_affine_in(log(x),{x}))
        ARIADNE_TEST_ASSERT(not is_affine_in(sqrt(x),{x}))
        ARIADNE_TEST_ASSERT(not is_affine_in(exp(x),{x}))
        ARIADNE_TEST_ASSERT(not is_affine_in(rec(x),{x}))
        ARIADNE_TEST_ASSERT(is_affine_in(hlf(x),{x}))
        ARIADNE_TEST_ASSERT(is_affine_in(nul(exp(x)),{x}))
        ARIADNE_TEST_ASSERT(not is_affine_in(sqr(x),{x}))
        ARIADNE_TEST_ASSERT(is_affine_in(neg(x),{x}))
        ARIADNE_TEST_ASSERT(not is_affine_in(pow(x,3),{x}));
        ARIADNE_TEST_ASSERT(not is_affine_in(sqr(x),{x}));
        ARIADNE_TEST_ASSERT(not is_affine_in(x*y,{x,y}));
        ARIADNE_TEST_ASSERT(not is_affine_in(x*x,{x}));
        ARIADNE_TEST_ASSERT(not is_affine_in(0*x*x,{x}));
        ARIADNE_TEST_ASSERT(not is_affine_in(x/y,{y}));
    }

    Void test_is_polynomial_in() {
        ARIADNE_TEST_ASSERT(is_polynomial_in(1,{x}))
        ARIADNE_TEST_ASSERT(is_polynomial_in(x,{x}))
        ARIADNE_TEST_ASSERT(is_polynomial_in(x*x,{x}))
        ARIADNE_TEST_ASSERT(not is_polynomial_in(atan(x),{x}))
        ARIADNE_TEST_ASSERT(not is_polynomial_in(acos(x),{x}))
        ARIADNE_TEST_ASSERT(not is_polynomial_in(asin(x),{x}))
        ARIADNE_TEST_ASSERT(not is_polynomial_in(cos(x),{x}))
        ARIADNE_TEST_ASSERT(not is_polynomial_in(sin(x),{x}))
        ARIADNE_TEST_ASSERT(not is_polynomial_in(tan(x),{x}))
        ARIADNE_TEST_ASSERT(not is_polynomial_in(abs(x),{x}))
        ARIADNE_TEST_ASSERT(not is_polynomial_in(log(x),{x}))
        ARIADNE_TEST_ASSERT(not is_polynomial_in(sqrt(x),{x}))
        ARIADNE_TEST_ASSERT(not is_polynomial_in(exp(x),{x}))
        ARIADNE_TEST_ASSERT(not is_polynomial_in(rec(x),{x}))
        ARIADNE_TEST_ASSERT(is_polynomial_in(hlf(x),{x}))
        ARIADNE_TEST_ASSERT(is_polynomial_in(nul(exp(x)),{x}))
        ARIADNE_TEST_ASSERT(is_polynomial_in(sqr(x),{x}))
        ARIADNE_TEST_ASSERT(is_polynomial_in(neg(x),{x}))
        ARIADNE_TEST_ASSERT(is_polynomial_in(pow(x,2),{x}))
        ARIADNE_TEST_ASSERT(is_polynomial_in(pow(x,3),{x}))
        ARIADNE_TEST_ASSERT(is_polynomial_in(pow(x+sqr(x),3),{x}))
        ARIADNE_TEST_ASSERT(is_polynomial_in(x/y,{x}))
        ARIADNE_TEST_ASSERT(not is_polynomial_in(x/y,{y}))
        ARIADNE_TEST_ASSERT(not is_polynomial_in(x/y,{x,y}))
        ARIADNE_TEST_ASSERT(is_polynomial_in({x,sqr(y)},{x,y}))
        ARIADNE_TEST_ASSERT(not is_polynomial_in({x/y,sqr(y)},{x,y}))
    }

    Void test_vector_expression()
    {
        ARIADNE_TEST_SAME_TYPE(RealVariables,Variables<Real>);
        ARIADNE_TEST_SAME_TYPE(RealVectorExpression,Expression<RealVector>);
        RealVector cs={-2,3,-5};
        RealVariables xs("x",3);
        RealVectorVariable vv("v",3);
        ARIADNE_TEST_CONSTRUCT(RealVectorExpression,vce,(cs));
        ARIADNE_TEST_CONSTRUCT(RealVectorExpression,vexs,(xs));
        ARIADNE_TEST_CONSTRUCT(RealVectorExpression,evv,(vv));
        ARIADNE_TEST_CONSTRUCT(RealVectorExpression,vle,({xs[0],cs[1],xs[2]+cs[2]}));
        ARIADNE_TEST_EQUALS(vce.size(),cs.size());
        ARIADNE_TEST_EQUALS(vexs.size(),xs.size());
        ARIADNE_TEST_EQUALS(evv.size(),3u);
        ARIADNE_TEST_EQUALS(vle.size(),3u);
        ARIADNE_TEST_EQUALS((vce+evv).size(),vce.size());
        RealVectorExpression vee=Vector<RealExpression>({xs[0],xs[1],xs[2]});
        RealVectorExpression vex=xs;
        ARIADNE_TEST_ASSERT(vex.node_ptr().operator->()==vex.node_raw_ptr());
        ARIADNE_TEST_ASSERT(vex.node_raw_ptr()==&vex.node_ref());
        RealVectorExpression vyce({cs[0],y,cs[2]});
        RealValuation val({xs[0]=7, xs[1]=11, xs[2]=13});
        RealVector vec({val[xs[0]], val[xs[1]], val[xs[2]]});
        RealValuation valv({vv[0]=7, vv[1]=11, vv[2]=13});
        RealVector vecv({valv[vv[0]], valv[vv[1]], valv[vv[2]]});
        ARIADNE_TEST_EQUAL(evaluate(vce,val),cs);
        ARIADNE_TEST_PRINT(evaluate(vex,val));
        ARIADNE_TEST_EQUAL(evaluate(vex,val),vec);
        ARIADNE_TEST_EQUAL(evaluate(vex+vce,val)[0],val[xs[0]]+cs[0]);
        ARIADNE_TEST_EQUAL(evaluate(vex+vce,val),evaluate(vex,val)+cs);
        ARIADNE_TEST_PRINT((vex+vce)/vex[1]);
        ARIADNE_TEST_PRINT(evaluate((vex+vce)/vex[1],val));
        ARIADNE_TEST_PRINT(simplify(((vex+vce)/vex[1])[2]));
        ARIADNE_TEST_SAME(simplify(((vex+vce)/vex[1])[2]),(xs[2]+cs[2])/xs[1]);

        ARIADNE_TEST_SAME(simplify((vex+vce)[0]),xs[0]+cs[0])
        ARIADNE_TEST_SAME(simplify((vex+vce)[1]),xs[1]+cs[1])
        ARIADNE_TEST_SAME(simplify((vex-vce)[0]),xs[0]-cs[0])
        ARIADNE_TEST_SAME(simplify((vex*y)[1]),xs[1]*y)
        ARIADNE_TEST_SAME(simplify((y*(vex-vce))[0]),y*(xs[0]-cs[0]))
        ARIADNE_TEST_SAME(simplify(((vex-vce)*y)[0]),(xs[0]-cs[0])*y)
        ARIADNE_TEST_SAME(simplify(((vex-vce)/y)[0]),(xs[0]-cs[0])/y)

    }

    Void test_function()
    {
        // Test to ensure that constants are handled correctly.
        Real tc=Dyadic(5.0);
        Real tx=Dyadic(1.125);
        Real ty=Dyadic(2.375);
        Real tz=Dyadic(3.750);

        Vector<Real> tv={tx,ty,tz};
        Valuation<Real> tw({{x,tx},{y,ty},{z,tz}});

        RealConstant c("5",tc);

        RealExpression e1=c;
        EffectiveScalarMultivariateFunction f1=make_function(RealSpace(List<RealVariable>({x,y,z})),e1);
        ARIADNE_TEST_PRINT(f1);
        ARIADNE_TEST_EQUAL(f1.evaluate(tv), tc);

        RealExpression e2=c+x;
        EffectiveScalarMultivariateFunction f2=make_function(RealSpace(List<RealVariable>({x,y,z})),e2);
        ARIADNE_TEST_PRINT(f2);
        ARIADNE_TEST_EQUAL(f2.evaluate(tv), tc+tx);

        RealExpression e3=c+x+c*y;
        EffectiveScalarMultivariateFunction f3=make_function({x,y,z},e3);
        ARIADNE_TEST_PRINT(f3);
        ARIADNE_TEST_EQUAL(f3.evaluate(tv), tc+tx+tc*ty);

        RealExpression e4=exp(c+x);
        EffectiveScalarMultivariateFunction f4=make_function({x,y,z},e4);
        ARIADNE_TEST_PRINT(f4);
        ARIADNE_TEST_ASSERT(possibly((f4.evaluate(tv) == exp(tc+tx)).check(Effort(0))));

        RealElementaryAlgebra ax=RealExpression(x);
        RealElementaryAlgebra ay=RealExpression(y);
        RealElementaryAlgebra az=RealExpression(z);
        Vector<RealElementaryAlgebra> va={ax,ay,az};
        ARIADNE_TEST_PRINT(va.element_characteristics());
        ARIADNE_TEST_PRINT(va.zero_element());
        ARIADNE_TEST_PRINT(va);
        ARIADNE_TEST_PRINT(f3(va));
        ARIADNE_TEST_CONSTRUCT(RealExpression,g3,(f3(va).extract<RealExpression>()));
        ARIADNE_TEST_PRINT(e3);
        ARIADNE_TEST_PRINT(f3);
        ARIADNE_TEST_PRINT(g3);
        ARIADNE_TEST_PRINT(5+va[0]+5*va[1]);
        ARIADNE_TEST_EQUALS(evaluate(g3,tw),evaluate(e3,tw));
    }

    Void test() {
        ARIADNE_TEST_CALL(test_variables());
        ARIADNE_TEST_CALL(test_space());
        ARIADNE_TEST_CALL(test_constants());
        ARIADNE_TEST_CALL(test_valuations());
        ARIADNE_TEST_CALL(test_templates());
        ARIADNE_TEST_CALL(test_scalar_expression_templates());
        ARIADNE_TEST_CALL(test_expression());
        ARIADNE_TEST_CALL(test_write());
        ARIADNE_TEST_CALL(test_assignment());
        ARIADNE_TEST_CALL(test_parameters());
        ARIADNE_TEST_CALL(test_identical());
        ARIADNE_TEST_CALL(test_derivative());
        ARIADNE_TEST_CALL(test_simplify());
        ARIADNE_TEST_CALL(test_ordering());
        ARIADNE_TEST_CALL(test_count_nodes());
        ARIADNE_TEST_CALL(test_count_distinct_nodes());
        ARIADNE_TEST_CALL(test_count_distinct_node_pointers());
        ARIADNE_TEST_CALL(test_eliminate_common_subexpressions());
        ARIADNE_TEST_CALL(test_substitute());
        ARIADNE_TEST_CALL(test_is_constant_in());
        ARIADNE_TEST_CALL(test_is_additive_in());
        ARIADNE_TEST_CALL(test_is_affine_in());
        ARIADNE_TEST_CALL(test_is_polynomial_in());
        ARIADNE_TEST_CALL(test_function());
        ARIADNE_TEST_CALL(test_vector_expression());
    }

};


Int main() {
    TestExpression().test();
    return ARIADNE_TEST_FAILURES;
}
