# Helpers every example's `reproduce_<name>(mode)` uses: typed config, output directory, metrics.

using Configurations
using TOML

const MODES = (:quick, :full)

"""
Parse `config/<name>.toml` as the `@option` type `T`, then `check_config` it. Parsing fails on an
unknown, missing or mistyped key, and both modes are checked whichever one runs, so a typo in
`[full]` is caught by a quick run.
"""
function load_config(::Type{T}, name::AbstractString) where {T}
    config = from_toml(T, joinpath(REPOSITORY_ROOT, "config", "$name.toml"))
    check_config(config)
    return config
end

"""`output/<mode>/<name>/`, created if missing."""
function output_dir(name::AbstractString, mode::Symbol)
    mode in MODES || throw(ArgumentError("mode must be :quick or :full"))
    return mkpath(joinpath(REPOSITORY_ROOT, "output", String(mode), name))
end

"""Throw an `ArgumentError` with `message` unless `condition` holds."""
require(condition::Bool, message::AbstractString) = condition || throw(ArgumentError(message))

"""The `metrics.toml` table name of an ε, such as `epsilon_0_02`."""
epsilon_key(epsilon) = "epsilon_" * replace(string(epsilon), "." => "_")

function save_toml(path::AbstractString, table::AbstractDict)
    open(path, "w") do io
        return TOML.print(io, table; sorted=true)
    end
    return path
end

"""
Write `metrics` to `<output>/metrics.toml`, headed by the config as run: `mode`, the `[model]`
table and that mode's table as `[settings]`, so every number can be traced to the settings
behind it.
"""
function save_metrics(output, mode::Symbol, model, settings, metrics::AbstractDict)
    config = Dict{String,Any}(
        "mode" => String(mode), "model" => to_dict(model), "settings" => to_dict(settings)
    )
    isempty(intersect(keys(config), keys(metrics))) ||
        throw(ArgumentError("metrics must not reuse the keys mode, model or settings"))
    return save_toml(joinpath(output, "metrics.toml"), merge(config, metrics))
end
