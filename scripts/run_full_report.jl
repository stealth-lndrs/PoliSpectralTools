#!/usr/bin/env julia

using Test
using Printf

const COUNT_KEYS = (:passes, :fails, :errors, :broken)
const CountTuple = NamedTuple{COUNT_KEYS, NTuple{4, Int}}

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))

function run_test_suite()
    return @testset "SpectralTools" begin
        include(joinpath(PROJECT_ROOT, "test", "runtests.jl"))
    end
end

function to_counts(tc::Test.TestCounts)::CountTuple
    return (passes = tc.passes, fails = tc.fails, errors = tc.errors, broken = tc.broken)
end

function sum_counts(a::CountTuple, b::CountTuple)::CountTuple
    return (passes = a.passes + b.passes,
            fails = a.fails + b.fails,
            errors = a.errors + b.errors,
            broken = a.broken + b.broken)
end

function aggregate_counts(ts::Test.DefaultTestSet)::CountTuple
    total = to_counts(Test.get_test_counts(ts))
    for result in ts.results
        result isa Test.DefaultTestSet || continue
        total = sum_counts(total, aggregate_counts(result))
    end
    return total
end

function collect_testsets(ts::Test.DefaultTestSet, successes, problematic; depth = 0)
    for result in ts.results
        result isa Test.DefaultTestSet || continue
        counts = aggregate_counts(result)
        entry = (repeat("  ", depth) * result.description, counts)
        if counts.fails == 0 && counts.errors == 0
            push!(successes, entry)
        else
            push!(problematic, entry)
        end
        collect_testsets(result, successes, problematic; depth = depth + 1)
    end
end

function collect_failures(ts::Test.DefaultTestSet, fails, errors, stack)
    for result in ts.results
        if result isa Test.DefaultTestSet
            collect_failures(result, fails, errors, vcat(stack, result.description))
        elseif result isa Test.Fail
            push!(fails, (copy(stack), result))
        elseif result isa Test.Error
            push!(errors, (copy(stack), result))
        end
    end
end

function show_result(io, label, counts::CountTuple)
    @printf(io, "  %-25s  pass=%4d  fail=%2d  error=%2d  broken=%2d\n",
            label, counts.passes, counts.fails, counts.errors, counts.broken)
end

function print_failure_block(io, entries, title)
    println(io, title)
    isempty(entries) && println(io, "  (none)\n")
    for (path, item) in entries
        location = join(path, " → ")
        println(io, "• ", isempty(location) ? "(root)" : location)
        src = getproperty(item, :source)
        if src !== nothing
            println(io, "    at $(src.file):$(src.line)")
        end
        buf = IOBuffer()
        show(buf, item)
        println(io, "    ", replace(String(take!(buf)), '\n' => "\n    "))
        println(io)
    end
end

function parse_cov_file(path)
    covered = 0
    total = 0
    open(path, "r") do io
        for line in eachline(io)
            token = strip(line)
            isempty(token) && continue
            token == "-" && continue
            total += 1
            val = maybe_tryparse_int(token)
            covered += (val !== nothing && val > 0) ? 1 : 0
        end
    end
    return covered, total
end

function maybe_tryparse_int(token)
    try
        return parse(Int, token)
    catch
        return nothing
    end
end

function coverage_summary()
    src_dir = joinpath(PROJECT_ROOT, "src")
    files = sort(filter(f -> endswith(f, ".jl"), readdir(src_dir; join = true)))
    entries = Tuple{String,Int,Int}[]
    for src in files
        cov_path = src * ".cov"
        isfile(cov_path) || continue
        covered, total = parse_cov_file(cov_path)
        push!(entries, (basename(src), covered, total))
    end
    return entries
end

function print_coverage(entries)
    println("Coverage Summary")
    if isempty(entries)
        println("  No coverage artifacts found.")
        println("  Run this script with `julia --project --code-coverage=user scripts/run_full_report.jl`.")
        return
    end
    total_cov = 0
    total_lines = 0
    for (name, covered, total) in entries
        total_lines += total
        total_cov += covered
        pct = total == 0 ? 100.0 : 100 * covered / total
        @printf("  %-20s %6.2f%% (%d/%d)\n", name, pct, covered, total)
    end
    pct_total = total_lines == 0 ? 100.0 : 100 * total_cov / total_lines
    @printf("  %-20s %6.2f%% (%d/%d)\n\n", "TOTAL", pct_total, total_cov, total_lines)
end

function main()
    ts = run_test_suite()
    counts = aggregate_counts(ts)
    status = (counts.fails == 0 && counts.errors == 0) ? "PASS" : "FAIL"
    println("=== SpectralTools – Full Package Test Report ===")
    println("Status: ", status)
    show_result(stdout, "Overall", counts)
    println()

    successes = Tuple{String,CountTuple}[]
    problematic = Tuple{String,CountTuple}[]
    collect_testsets(ts, successes, problematic)
    println("Testsets passing cleanly:")
    isempty(successes) && println("  (none)")
    for (label, ct) in successes
        show_result(stdout, label, ct)
    end
    println()
    println("Testsets with issues:")
    isempty(problematic) && println("  (none)")
    for (label, ct) in problematic
        show_result(stdout, label, ct)
    end
    println()

    fails = Vector{Tuple{Vector{String},Test.Fail}}()
    errors = Vector{Tuple{Vector{String},Test.Error}}()
    collect_failures(ts, fails, errors, [ts.description])
    print_failure_block(stdout, fails, "Failures:")
    print_failure_block(stdout, errors, "Errors:")

    entries = coverage_summary()
    println()
    print_coverage(entries)
end

main()
