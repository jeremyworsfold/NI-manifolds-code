#!/usr/bin/env julia

include(joinpath(@__DIR__, "lib", "common.jl"))

# Every example feeds a manuscript figure.
const EXAMPLES = ("tk2", "ants")

function parse_arguments(arguments)
    mode = :quick
    example = "all"
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        if argument == "--mode"
            index == length(arguments) && error("--mode requires quick or full")
            mode = Symbol(arguments[index + 1])
            index += 2
        elseif argument == "--example"
            index == length(arguments) && error("--example requires a model name or all")
            example = arguments[index + 1]
            index += 2
        else
            error("unknown argument: $argument")
        end
    end
    mode in (:quick, :full) || error("--mode must be quick or full")
    example in ("all", EXAMPLES...) ||
        error("--example must be one of: all, $(join(EXAMPLES, ", "))")
    return mode, example
end

function main(arguments=ARGS)
    mode, example = parse_arguments(arguments)
    selected = example == "all" ? EXAMPLES : (example,)
    for item in selected
        # Each script `<item>.jl` defines `reproduce_<item>(mode)`.
        include(joinpath(@__DIR__, "$item.jl"))
        reproduce = Base.invokelatest(getglobal, Main, Symbol("reproduce_$item"))
        output = Base.invokelatest(reproduce, mode)
        println("generated $item artifacts in $output")
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
