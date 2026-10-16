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

#include "drawing/embedded.hpp"
#ifdef ARIADNE_HAS_CAIRO
#include "drawing/cairo.hpp"
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
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
    if(text.empty() || text.front()=='-') throw std::invalid_argument("invalid positive integer: "+text);
    std::size_t parsed=0;
    const auto value=std::stoull(text,&parsed,10);
    if(parsed!=text.size() || value==0 || value>limit)
        throw std::invalid_argument("integer out of range: "+text);
    return static_cast<unsigned>(value);
}

void usage() {
    std::cout <<
        "Usage: benchmark_drawing [--backend embedded|cairo|both] [--scene NAME|all]\n"
        "       [--width N] [--height N] [--iterations N] [--warmup N]\n"
        "       [--write-png] [--png-prefix PATH]\n"
        "Scenes: rectangles, polylines, circles, polygons, winding, mixed\n"
        "Outputs CSV: absolute render timings exclude canvas construction,\n"
        "initialisation (including axis labels), and optional PNG I/O.\n";
}

Options parse(int argc,char const* const* argv) {
    Options result;
    for(int i=1;i<argc;++i) {
        const std::string arg(argv[i]);
        if(arg=="--help" || arg=="-h") { usage(); std::exit(0); }
        if(arg=="--write-png") { result.write_png=true; continue; }
        if(i+1>=argc) throw std::invalid_argument("missing value for "+arg);
        const std::string value(argv[++i]);
        if(arg=="--backend") result.backend=value;
        else if(arg=="--scene") result.scene=value;
        else if(arg=="--width") result.width=positive_integer(value,16384);
        else if(arg=="--height") result.height=positive_integer(value,16384);
        else if(arg=="--iterations") result.iterations=positive_integer(value,10000);
        else if(arg=="--warmup") result.warmup=positive_integer(value,10000);
        else if(arg=="--png-prefix") result.png_prefix=value;
        else throw std::invalid_argument("unknown argument: "+arg);
    }
    if(result.backend!="embedded" && result.backend!="cairo" && result.backend!="both")
        throw std::invalid_argument("invalid backend: "+result.backend);
#ifndef ARIADNE_HAS_CAIRO
    if(result.backend!="embedded")
        throw std::invalid_argument("Cairo was not compiled; use --backend embedded");
#endif
    if(result.scene!="all" && result.scene!="rectangles" &&
       result.scene!="polylines" && result.scene!="circles" &&
       result.scene!="polygons" && result.scene!="winding" &&
       result.scene!="mixed")
        throw std::invalid_argument("invalid scene: "+result.scene);
    // Avoid overflowing the underlying image surface allocation.
    if(static_cast<std::uint64_t>(result.width+170u)*(result.height+50u)>100000000ULL)
        throw std::invalid_argument("canvas too large");
    return result;
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

void benchmark(Options const& opts,std::string const& backend,Scene const& scene) {
    Samples samples;
    const unsigned repeats=opts.warmup+opts.iterations;
    const std::string filename=opts.png_prefix+"-"+backend+"-"+scene.name+".png";
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
    }
    const auto c=summarize(samples.creation);
    const auto init=summarize(samples.initialisation);
    const auto render=summarize(samples.render);
    const auto write=summarize(samples.export_png);
    const double throughput=render.median>0.0 ?
        1000.0*static_cast<double>(scene.primitives)/render.median : 0.0;
    std::cout << backend << ',' << scene.name << ','
              << opts.width << ',' << opts.height << ','
              << scene.primitives << ',' << scene.commands.size() << ','
              << opts.iterations << ','
              << std::fixed << std::setprecision(4)
              << c.median << ',' << init.median << ','
              << render.median << ',' << render.p95 << ','
              << throughput << ',';
    if(opts.write_png) std::cout << write.median;
    std::cout << '\n';
}

} // namespace

int main(int argc,char const* const* argv) {
    try {
        const auto opts=parse(argc,argv);
        const auto scenes=make_scenes(); // Scene geometry is precomputed, never timed.
        std::cout << "backend,scene,width,height,primitives,commands,iterations,"
                     "create_median_ms,initialise_median_ms,render_median_ms,"
                     "render_p95_ms,primitives_per_second,png_write_median_ms\n";
        for(auto const& scene:scenes) {
            if(opts.scene!="all" && scene.name!=opts.scene) continue;
            if(opts.backend=="embedded" || opts.backend=="both")
                benchmark(opts,"embedded",scene);
            if(opts.backend=="cairo" || opts.backend=="both")
                benchmark(opts,"cairo",scene);
        }
        return 0;
    } catch(std::exception const& e) {
        std::cerr << "benchmark_drawing: " << e.what() << '\n';
        usage();
        return 1;
    }
}
