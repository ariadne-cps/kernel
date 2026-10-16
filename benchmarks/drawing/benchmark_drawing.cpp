/***************************************************************************
 *            benchmarks/drawing/benchmark_drawing.cpp
 *
 *  Copyright 2026 Luca Geretti
 *
 *  This file is part of Ariadne, released under the GPL v3 or later.
 ***************************************************************************/

// This benchmark deliberately drives CanvasInterface directly. In particular,
// no Figure, Drawer, constrained set, approximation, or geometry computation
// is included in the measured rendering time.

#include "drawing/command_line_interface.hpp"
#include "drawing/embedded.hpp"
#include "logging/logging.hpp"
#ifdef ARIADNE_HAS_CAIRO
#include "drawing/cairo.hpp"
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Ariadne;

namespace {

using Clock = std::chrono::steady_clock;
constexpr double pi = 3.14159265358979323846;

enum class Operation {
    move, line, circle, fill, stroke, width, line_colour, fill_colour, opacity
};

struct Command {
    Operation op;
    double a=0.0, b=0.0, c=0.0;
};

struct Scene {
    std::string name;
    std::vector<Command> commands;
    std::size_t primitives=0;

    void add(Operation op,double a=0.0,double b=0.0,double c=0.0) {
        commands.push_back({op,a,b,c});
    }
    void rectangle(double xl,double yl,double xu,double yu) {
        add(Operation::move,xl,yl);
        add(Operation::line,xu,yl);
        add(Operation::line,xu,yu);
        add(Operation::line,xl,yu);
        add(Operation::line,xl,yl);
        add(Operation::fill);
        ++primitives;
    }
    void circle(double x,double y,double radius) {
        add(Operation::circle,x,y,radius);
        add(Operation::fill);
        ++primitives;
    }
};

struct Options {
    unsigned width=640, height=480;
    unsigned iterations=5, warmup=1;
    std::string backend="embedded", scene="all", png_prefix="benchmark-drawing";
    bool write_png=false;
};

unsigned positive_integer(std::string const& text,unsigned limit) {
    if(text.empty() || text.front()=='-')
        throw std::invalid_argument("invalid positive integer: "+text);
    std::size_t parsed=0;
    const auto value=std::stoull(text,&parsed,10);
    if(parsed!=text.size() || value==0 || value>limit)
        throw std::invalid_argument("integer out of range: "+text);
    return static_cast<unsigned>(value);
}

List<ArgumentParser> benchmark_argument_parsers(Options& opts) {
    List<ArgumentParser> parsers;
    parsers.append(make_valued_argument_parser(
        "", "backend", "Drawing backend: embedded, cairo, or both (default: embedded)",
        [&opts](String const& input) -> VoidFunction {
            const std::string value=input;
            if(value!="embedded" && value!="cairo" && value!="both")
                throw std::invalid_argument("invalid drawing backend: "+value);
#ifndef ARIADNE_HAS_CAIRO
            if(value!="embedded")
                throw std::invalid_argument("optional reference backend not compiled");
#endif
            return [&opts,value] { opts.backend=value; };
        }));
    parsers.append(make_valued_argument_parser(
        "", "scene", "Scenario: rectangles, polylines, circles, polygons, winding, mixed, or all",
        [&opts](String const& input) -> VoidFunction {
            const std::string value=input;
            if(value!="all" && value!="rectangles" && value!="polylines" &&
               value!="circles" && value!="polygons" && value!="winding" &&
               value!="mixed")
                throw std::invalid_argument("invalid drawing scene: "+value);
            return [&opts,value] { opts.scene=value; };
        }));
    parsers.append(make_valued_argument_parser(
        "", "width", "Plot width in pixels (1..16384; default: 640)",
        [&opts](String const& input) -> VoidFunction {
            const auto value=positive_integer(input,16384u);
            return [&opts,value] { opts.width=value; };
        }));
    parsers.append(make_valued_argument_parser(
        "", "height", "Plot height in pixels (1..16384; default: 480)",
        [&opts](String const& input) -> VoidFunction {
            const auto value=positive_integer(input,16384u);
            return [&opts,value] { opts.height=value; };
        }));
    parsers.append(make_valued_argument_parser(
        "", "iterations", "Measured repetitions per scenario (1..10000; default: 5)",
        [&opts](String const& input) -> VoidFunction {
            const auto value=positive_integer(input,10000u);
            return [&opts,value] { opts.iterations=value; };
        }));
    parsers.append(make_valued_argument_parser(
        "", "warmup", "Warmup repetitions per scenario (1..10000; default: 1)",
        [&opts](String const& input) -> VoidFunction {
            const auto value=positive_integer(input,10000u);
            return [&opts,value] { opts.warmup=value; };
        }));
    parsers.append(make_unvalued_argument_parser(
        "", "write-png", "Measure optional PNG encoding and filesystem output",
        [&opts] { opts.write_png=true; }));
    parsers.append(make_valued_argument_parser(
        "", "png-prefix", "Output path prefix for optional PNG files",
        [&opts](String const& input) -> VoidFunction {
            const std::string value=input;
            if(value.empty()) throw std::invalid_argument("empty PNG prefix");
            return [&opts,value] { opts.png_prefix=value; };
        }));
    return parsers;
}

void validate_options(Options const& opts) {
    // Include reserved plot margins when checking memory limits.
    if(static_cast<std::uint64_t>(opts.width+170u)*(opts.height+50u)>100000000ULL)
        throw std::invalid_argument("canvas too large");
}

void set_style(Scene& scene,double width,double opacity) {
    scene.add(Operation::width,width);
    scene.add(Operation::line_colour,0.08,0.12,0.36);
    scene.add(Operation::fill_colour,0.76,0.34,0.16);
    scene.add(Operation::opacity,opacity);
}

Scene rectangles() {
    Scene scene{"rectangles",{},0};
    set_style(scene,0.8,0.6);
    constexpr int columns=24,rows=16;
    for(int iy=0;iy<rows;++iy) for(int ix=0;ix<columns;++ix) {
        const double dx=1.8/columns, dy=1.8/rows;
        const double x=-0.9+ix*dx, y=-0.9+iy*dy;
        scene.rectangle(x,y,x+dx*0.82,y+dy*0.82);
    }
    return scene;
}

Scene polylines() {
    Scene scene{"polylines",{},0};
    set_style(scene,2.0,1.0);
    constexpr int paths=72, segments=24;
    for(int j=0;j<paths;++j) {
        const double offset=-0.85+1.7*j/(paths-1);
        scene.add(Operation::move,-0.95,offset);
        for(int k=1;k<=segments;++k) {
            const double x=-0.95+1.9*k/segments;
            const double y=offset+0.025*std::sin(7.0*x+0.19*j);
            scene.add(Operation::line,x,y);
        }
        scene.add(Operation::stroke);
        ++scene.primitives;
    }
    return scene;
}

Scene circles() {
    Scene scene{"circles",{},0};
    set_style(scene,1.2,0.55);
    constexpr int columns=12,rows=10;
    for(int iy=0;iy<rows;++iy) for(int ix=0;ix<columns;++ix) {
        const double x=-0.85+1.7*ix/(columns-1);
        const double y=-0.83+1.66*iy/(rows-1);
        const double radius=0.032+0.004*((ix+iy)%3);
        scene.circle(x,y,radius);
    }
    return scene;
}

Scene polygons() {
    Scene scene{"polygons",{},0};
    set_style(scene,1.0,0.65);
    constexpr int columns=10,rows=8,vertices=12;
    for(int iy=0;iy<rows;++iy) for(int ix=0;ix<columns;++ix) {
        const double cx=-0.82+1.64*ix/(columns-1);
        const double cy=-0.80+1.60*iy/(rows-1);
        for(int k=0;k<=vertices;++k) {
            const double angle=2*pi*(k%vertices)/vertices;
            const double radius=(k%2==0) ? 0.075 : 0.038;
            const double x=cx+radius*std::cos(angle);
            const double y=cy+radius*std::sin(angle);
            scene.add(k==0 ? Operation::move : Operation::line,x,y);
        }
        scene.add(Operation::fill);
        ++scene.primitives;
    }
    return scene;
}

Scene winding() {
    Scene scene{"winding",{},0};
    set_style(scene,1.0,0.7);
    constexpr int columns=9,rows=6;
    for(int iy=0;iy<rows;++iy) for(int ix=0;ix<columns;++ix) {
        const double x=-0.85+1.7*ix/(columns-1);
        const double y=-0.82+1.64*iy/(rows-1);
        const double outer=0.071,inner=0.037;
        // Outer CCW and inner CW in the same fill: requires nonzero winding.
        scene.add(Operation::move,x-outer,y-outer);
        scene.add(Operation::line,x+outer,y-outer);
        scene.add(Operation::line,x+outer,y+outer);
        scene.add(Operation::line,x-outer,y+outer);
        scene.add(Operation::line,x-outer,y-outer);
        scene.add(Operation::move,x-inner,y-inner);
        scene.add(Operation::line,x-inner,y+inner);
        scene.add(Operation::line,x+inner,y+inner);
        scene.add(Operation::line,x+inner,y-inner);
        scene.add(Operation::line,x-inner,y-inner);
        scene.add(Operation::fill);
        ++scene.primitives;
    }
    return scene;
}

std::vector<Scene> make_scenes() {
    std::vector<Scene> scenes;
    scenes.push_back(rectangles());
    scenes.push_back(polylines());
    scenes.push_back(circles());
    scenes.push_back(polygons());
    scenes.push_back(winding());
    Scene mixed{"mixed",{},0};
    for(auto const& scene:scenes) {
        mixed.commands.insert(mixed.commands.end(),scene.commands.begin(),scene.commands.end());
        mixed.primitives+=scene.primitives;
    }
    scenes.push_back(std::move(mixed));
    return scenes;
}

void draw(CanvasInterface& canvas,Scene const& scene) {
    for(auto const& command:scene.commands) {
        switch(command.op) {
            case Operation::move: canvas.move_to(command.a,command.b); break;
            case Operation::line: canvas.line_to(command.a,command.b); break;
            case Operation::circle: canvas.circle(command.a,command.b,command.c); break;
            case Operation::fill: canvas.fill(); break;
            case Operation::stroke: canvas.stroke(); break;
            case Operation::width: canvas.set_line_width(command.a); break;
            case Operation::line_colour:
                canvas.set_line_colour(command.a,command.b,command.c); break;
            case Operation::fill_colour:
                canvas.set_fill_colour(command.a,command.b,command.c); break;
            case Operation::opacity: canvas.set_fill_opacity(command.a); break;
            default: throw std::logic_error("invalid drawing benchmark command");
        }
    }
}

SharedPointer<CanvasInterface> make_canvas(std::string const& backend,
                                           unsigned width,unsigned height) {
    if(backend=="embedded")
        return EmbeddedGraphicsBackend().make_canvas("",width,height,false);
#ifdef ARIADNE_HAS_CAIRO
    if(backend=="cairo")
        return CairoGraphicsBackend().make_canvas("",width,height,false);
#endif
    throw std::invalid_argument("unsupported backend: "+backend);
}

double milliseconds(Clock::time_point a,Clock::time_point b) {
    return std::chrono::duration<double,std::milli>(b-a).count();
}

struct Samples {
    std::vector<double> creation,initialisation,render,export_png;
};

struct Statistics {
    double median=0.0,p95=0.0;
};
Statistics summarize(std::vector<double> values) {
    if(values.empty()) return {};
    std::sort(values.begin(),values.end());
    const std::size_t middle=values.size()/2;
    const double median=values.size()%2 ? values[middle] :
        (values[middle-1]+values[middle])/2.0;
    const std::size_t p95_index=static_cast<std::size_t>(
        std::ceil(0.95*static_cast<double>(values.size())))-1u;
    return {median,values[p95_index]};
}

struct Result {
    std::string backend;
    std::string scene;
    std::size_t primitives=0, commands=0;
    Statistics creation, initialisation, rendering, png;
    double primitives_per_second=0.0;
};

Result benchmark(Options const& opts,std::string const& backend,
                 Scene const& scene) {
    Samples samples;
    const unsigned repeats=opts.warmup+opts.iterations;
    const std::string filename=
        opts.png_prefix+"-"+backend+"-"+scene.name+".png";

    LOGGING_PRINTLN_AT(2, "Measuring backend="<<backend
        <<", scene="<<scene.name<<", iterations="<<opts.iterations)

    for(unsigned i=0;i<repeats;++i) {
        const auto before=Clock::now();
        auto canvas=make_canvas(backend,opts.width,opts.height);
        const auto constructed=Clock::now();
        canvas->initialise("x","y",-1.0,1.0,-1.0,1.0);
        const auto initialized=Clock::now();
        draw(*canvas,scene);
        const auto rendered=Clock::now();
        if(opts.write_png) canvas->write(filename.c_str());
        const auto finished=Clock::now();
        if(i>=opts.warmup) {
            samples.creation.push_back(milliseconds(before,constructed));
            samples.initialisation.push_back(milliseconds(constructed,initialized));
            samples.render.push_back(milliseconds(initialized,rendered));
            if(opts.write_png)
                samples.export_png.push_back(milliseconds(rendered,finished));
        }
        LOGGING_PRINTLN_AT(3, "  "<<backend<<"/"<<scene.name
            <<" sample="<<i+1u<<"/"<<repeats
            <<", render="<<milliseconds(initialized,rendered)<<" ms")
    }

    const auto creation=summarize(samples.creation);
    const auto initialisation=summarize(samples.initialisation);
    const auto rendering=summarize(samples.render);
    const auto png=summarize(samples.export_png);
    const double throughput=rendering.median>0.0 ?
        1000.0*static_cast<double>(scene.primitives)/rendering.median : 0.0;
    return {backend,scene.name,scene.primitives,scene.commands.size(),
            creation,initialisation,rendering,png,throughput};
}

// Consistent terminal formatting shared by all scenarios; no CSV is emitted.
// The overall report is built after timing to avoid logging in hot loops.
std::string fixed(double value,int precision=3) {
    std::ostringstream stream;
    stream<<std::fixed<<std::setprecision(precision)<<value;
    return stream.str();
}

void print_table(Options const& opts,std::vector<Result> const& results) {
    struct Column { std::string title; std::size_t width; bool left; };
    std::vector<Column> columns{
        {"backend",7u,true},
        {"scene",10u,true},
        {"shapes",6u,false},
        {"commands",8u,false},
        {"create ms",9u,false},
        {"init ms",8u,false},
        {"render ms",9u,false},
        {"p95 ms",8u,false},
        {"shapes/s",9u,false}
    };
    if(opts.write_png) columns.push_back({"PNG ms",8u,false});

    std::vector<std::vector<std::string>> rows;
    rows.reserve(results.size());
    for(auto const& result:results) {
        std::vector<std::string> row{
            result.backend,
            result.scene,
            std::to_string(result.primitives),
            std::to_string(result.commands),
            fixed(result.creation.median),
            fixed(result.initialisation.median),
            fixed(result.rendering.median),
            fixed(result.rendering.p95),
            fixed(result.primitives_per_second,0)
        };
        if(opts.write_png) row.push_back(fixed(result.png.median));
        for(std::size_t i=0;i<row.size();++i)
            columns[i].width=std::max(columns[i].width,row[i].size());
        rows.push_back(std::move(row));
    }

    const std::size_t spacing=2u;
    std::size_t total_width=spacing*(columns.size()-1u);
    for(auto const& column:columns)
        total_width+=column.width;
    const std::string line(total_width,'-');

    std::ostringstream out;
    out<<"Drawing backend benchmark\n"
       <<"Plot: "<<opts.width<<" x "<<opts.height
       <<" px | iterations: "<<opts.iterations
       <<" | warmup: "<<opts.warmup
       <<" | scene: "<<opts.scene<<"\n"
       <<"Render excludes canvas creation, initialization and PNG I/O.\n\n";
    auto print_row=[&](std::vector<std::string> const& values) {
        for(std::size_t i=0;i<columns.size();++i) {
            if(i!=0) out<<std::string(spacing,' ');
            out<<(columns[i].left ? std::left : std::right)
               <<std::setw(static_cast<int>(columns[i].width))
               <<values[i];
        }
        out<<"\n";
    };
    std::vector<std::string> titles;
    titles.reserve(columns.size());
    for(auto const& column:columns) titles.push_back(column.title);
    print_row(titles);
    out<<line<<"\n";
    for(std::size_t i=0;i<rows.size();++i) {
        print_row(rows[i]);
        if(i+1u==rows.size() || rows[i+1u][1]!=rows[i][1])
            out<<line<<"\n";
    }
    LOGGING_PRINTLN(out.str())
}

} // namespace

int main(int argc,char const* const* argv) {
    Logger::instance().configuration().set_verbosity(1u);

    Options opts;
    auto command_line=make_kernel_command_line_interface();
    command_line.add(CommandLineContribution(
        "drawing-benchmark",benchmark_argument_parsers(opts)));
    if(not command_line.acquire(argc,argv))
        return -1;

    try {
        validate_options(opts);
        // All input geometry is prepared outside the measured rendering.
        const auto scenes=make_scenes();
        std::vector<Result> results;
        results.reserve(scenes.size()*(opts.backend=="both" ? 2u : 1u));
        for(auto const& scene:scenes) {
            if(opts.scene!="all" && scene.name!=opts.scene) continue;
            if(opts.backend=="embedded" || opts.backend=="both")
                results.push_back(benchmark(opts,"embedded",scene));
            if(opts.backend=="cairo" || opts.backend=="both")
                results.push_back(benchmark(opts,"cairo",scene));
        }
        print_table(opts,results);
        return 0;
    } catch(std::exception const& error) {
        LOGGING_PRINTLN("drawing benchmark failed: "<<error.what())
        return 1;
    }
}
