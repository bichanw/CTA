import PPSeq
const seq = PPSeq

# Other Imports
import PyPlot: plt
import DelimitedFiles: readdlm
import Random
import StatsBase: quantile

println("Runnign CTA/test_code.jl")
flush(stdout)


# read spikes
using MAT: matopen
file = matopen("/usr/people/bichanw/SpikeSorting/Codes/CTA/mat/julia.mat")
clu = read(file,"clu") # note that this does NOT introduce a variable ``varname`` into scope
t   = read(file,"t") # note that this does NOT introduce a variable ``varname`` into scope
filename = read(file,"filename")
close(file)
spikes = seq.Spike[]
for i in 1:length(clu)
    push!(spikes, seq.Spike(clu[i], t[i]))
end
# metadata
num_neurons = length(unique(clu));
max_time = t[end];


# plot raw raster
fig = seq.plot_raster(spikes; color=[0,0,0,0.05]) # returns matplotlib Figure
fig.set_size_inches([7, 3]);
fig.savefig("results/origianl_spikes.png")



config = Dict(

    # Model hyperparameters
    :num_sequence_types =>  2,
    :seq_type_conc_param => 1.0,
    :seq_event_rate => 0.04,

    :mean_event_amplitude => 100.0,
    :var_event_amplitude => 1000.0,
    
    :neuron_response_conc_param => 0.1,
    :neuron_offset_pseudo_obs => 1.0,
    :neuron_width_pseudo_obs => 1.0,
    :neuron_width_prior => 0.5,
    
    :num_warp_values => 1,
    :max_warp => 1.0,
    :warp_variance => 1.0,

    :mean_bkgd_spike_rate => 150.0, # is this for all neurons
    :var_bkgd_spike_rate => 50.0,
    :bkgd_spikes_conc_param => 0.3,
    :max_sequence_length => Inf,
    
    # # MCMC Sampling parameters.
    # :num_anneals => 10,
    # :samples_per_anneal => 100,
    # :max_temperature => 40.0,
    # :save_every_during_anneal => 10,
    # :samples_after_anneal => 2000,
    # :save_every_after_anneal => 10,
    # :split_merge_moves_during_anneal => 10,
    # :split_merge_moves_after_anneal => 10,
    # :split_merge_window => 1.0,

    # distributed MCMC
    :num_threads => 10,  # <--- This is the key parameter to include if you want to run parallel MCMC
    :num_anneals => 10,
    :samples_per_anneal => 100,
    :max_temperature => 40.0,
    :save_every_during_anneal => 10,
    :samples_after_anneal => 2000,
    :save_every_after_anneal => 10,
    :split_merge_moves_during_anneal => 0,  # SPLIT / MERGE not implemented for distributed MCMC
    :split_merge_moves_after_anneal => 0,   # SPLIT / MERGE not implemented for distributed MCMC
    :split_merge_window => 1.0,

);

# Initialize all spikes to background process.
init_assignments = fill(-1, length(spikes))

# Construct model struct (PPSeq instance).
model = seq.construct_model(config, max_time, num_neurons)

# Run Gibbs sampling with an initial annealing period.
results = seq.easy_sample!(model, spikes, init_assignments, config);

# Grab the final MCMC sample
final_globals = results[:globals_hist][end]
final_events = results[:latent_event_hist][end]
final_assignments = results[:assignment_hist][:, end]



# save data in matrix (matlab does not have PPSeq sturcture)
using JLD2
pwd()
jldsave("mat/"*filename*"_assignments.jld2";final_assignments,final_events)


# Helpful utility function that sorts the neurons to reveal sequences.
neuron_ordering = seq.sortperm_neurons(final_globals)

# Plot model-annotated raster.
fig = seq.plot_raster(
    spikes,
    final_events,
    final_assignments,
    neuron_ordering;
    color_cycle=["red", "blue"] # colors for each sequence type can be modified.
)
fig.set_size_inches([7, 3]);
fig.savefig("results/colored_spikes.png")

seq.plot_log_likes(config, results);

seq.plot_num_seq_events(config, results);

# # Create discrete time grid.
# num_timebins = 1000
# dt = max_time / num_timebins
# timebins = collect((0.5 * dt):dt:max_time)

# Compute a matrix firing rates (num_neurons x num_timebins)
F = seq.firing_rates(
    final_globals,
    final_events,
    timebins
)

# # Plot firing rates as a heatmap
# plt.imshow(F[neuron_ordering, :]; aspect="auto", origin="lower")
# plt.title("Firing Rates, last MCMC sample (spikes / second)")
# plt.ylabel("Neurons")
# plt.xlabel("timebins")
# plt.colorbar();

F_nrm = copy(F)
for n in 1:num_neurons
    F_nrm[n, :] .-= minimum(F[n, :])
    F_nrm[n, :] ./= maximum(F[n, :])
end

plt.title("Firing Rates, last MCMC sample (normalized)")
plt.ylabel("Neurons")
plt.xlabel("timebins")
plt.imshow(F_nrm[neuron_ordering, :]; aspect="auto", origin="lower")
plt.colorbar();
plt.savefig("results/FR.png");

# F_avg = zeros(num_neurons, num_timebins)

# for (G, E) in zip(results[:globals_hist], results[:latent_event_hist])
#     F_avg += seq.firing_rates(G, E, timebins)
# end
# F_avg ./= length(results[:globals_hist])

# # Plot average firing rates as a heatmap
# plt.imshow(F_avg[neuron_ordering, :]; aspect="auto", origin="lower")
# plt.title("Posterior Expected Firing Rates (spikes / second)")
# plt.ylabel("Neurons")
# plt.xlabel("timebins")
# plt.colorbar();

# F_avg_nrm = copy(F_avg)
# for n in 1:num_neurons
#     F_avg_nrm[n, :] .-= minimum(F_avg[n, :])
#     F_avg_nrm[n, :] ./= maximum(F_avg[n, :])
# end

# # Plot average firing rates as a heatmap
# plt.imshow(F_avg_nrm[neuron_ordering, :]; aspect="auto", origin="lower")
# plt.title("Posterior Expected Firing Rates (normalized)")
# plt.ylabel("Neurons")
# plt.xlabel("timebins")
# plt.colorbar();


