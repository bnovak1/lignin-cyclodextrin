from pathlib import Path
import numpy as np


configfile: "plumed.json"


output_dir = Path("../output/{lignol}/one_BCD")
analysis_dir = Path("../analysis/{lignol}/one_BCD")

# Copy some figures to the writing directory.
writing_dir = Path(
    "..", 
    "writing","HDBSCAN_clustering_for_identification_of_states_with_lignin_dimers_bound_to_cyclodextrin_from_molecular_dynamics_simulations",
    "000_Attachments"
)

# Need to get masses into pdb occupancy column for PLUMED.
# Note that masses guessed by MDAnalysis do not match masses from topology exactly.
rule nowater_pdb:
    input:
        gro=output_dir / "nowater.gro",
    output:
        pdb=str(Path("../analysis/{lignol}/nowater.pdb")),
    run:
        import MDAnalysis as mda

        u = mda.Universe(input.gro, input.gro)
        u.add_TopologyAttr("occupancies")
        u.atoms.occupancies = u.atoms.masses
        u.atoms.write(output.pdb)


rule nowater_pdbs:
    input:
        expand(
            rules.nowater_pdb.output.pdb,
            lignol=config["LIGNOLS"],
        ),


# Compute collective variables using PLUMED.
rule plumed_CVs:
    input:
        plumed=str(Path("{lignol}/CVs.plumed.dat")),
        traj=str(output_dir / "nowater.xtc"),
        pdb=rules.nowater_pdb.output.pdb,
    output:
        colvar=str(analysis_dir / "colvar.dat"),
    shell:
        "singularity exec ~/singularity/GROMACS_2018.3-PLUMED_2.6/ plumed driver --plumed {input.plumed} --timestep 4 --mf_xtc {input.traj} --pdb {input.pdb}"


rule plumed_CVss:
    input:
        expand(rules.plumed_CVs.output, lignol=config["LIGNOLS"]),


# Plot 2D and 3D scatter plots of the collective variables.
rule plots_2D:
    input:
        colvar=rules.plumed_CVs.output.colvar,
        script=Path("../scripts/plots.py"),
    output:
        scatter_dnorm_dtang=str(analysis_dir / "dnorm_dtang.png"),
        kde_dnorm_dtang=str(analysis_dir / "dnorm_dtang_kde.png"),
        scatter_dnorm_orient=str(analysis_dir / "dnorm_orient.png"),
        scatter_dtang_orient=str(analysis_dir / "dtang_orient.png"),
    params:
        output_dir=str(analysis_dir),
    shell:
        "python {input.script} --lignol {wildcards.lignol} --colvar {input.colvar} --outdir {params.output_dir} --dimension 2"


rule plots_2Ds:
    input:
        expand(rules.plots_2D.output, lignol=config["LIGNOLS"]),


rule plots_3D:
    input:
        colvar=rules.plumed_CVs.output.colvar,
        script=Path("../scripts/plots.py"),
    output:
        scatter_html=str(analysis_dir / "high_density_scatter.html"),
        scatter_png=str(analysis_dir / "high_density_scatter.png"),
    params:
        output_dir=str(analysis_dir),
    shell:
        """
        python {input.script} \
            --lignol {wildcards.lignol} \
            --colvar {input.colvar} \
            --outdir {params.output_dir} \
            --dimension 3
        convert {output.scatter_png} -trim {output.scatter_png}
        """


rule plots_3Ds:
    input:
        expand(rules.plots_3D.output, lignol=config["LIGNOLS"]),


# rule elbow_plot:
#     input:
#         colvar = rules.plumed_CVs.output.colvar,
#         script = Path("../scripts/elbow.py")
#     output:
#         elbow = str(analysis_dir / "elbow.png"),
#         max_eps = str(analysis_dir / "max_eps.txt")
#     params:
#         output_dir = str(analysis_dir)
#     shell:
#         "python {input.script} --lignol {wildcards.lignol} --colvar {input.colvar} --outdir {params.output_dir}"

# rule elbow_plots:
#     input:
#         expand(rules.elbow_plot.output, lignol=config['LIGNOLS'])

# min_samples_range = range(config["MIN_SAMPLES_RANGE"][0], config["MIN_SAMPLES_RANGE"][1]+1)
# min_samples_range_GG_BB = range(config["MIN_SAMPLES_RANGE_GG_BB"][0], config["MIN_SAMPLES_RANGE_GG_BB"][1]+1)

# HDBSCAN clustering.
rule clustering:
    input:
        colvar=rules.plumed_CVs.output.colvar,
        script=Path("../scripts/clustering.py"),
        read_script=Path("../scripts/read_colvar_file.py"),
    output:
        nclust=str(
            analysis_dir / "clustering/HDBSCAN_model/nclusters_{min_samples}.dat"
        ),
    threads: 8
    shell:
        "python {input.script} --colvar {input.colvar} --min_samples {wildcards.min_samples} --nclusters_file {output.nclust}"

rule clusterings:
    input:
        [
            f"../analysis/{lignol}/one_BCD/clustering/HDBSCAN_model/nclusters_{min_samples}.dat"
            for lignol in config["LIGNOLS"]
            for min_samples in range(
                config["MIN_SAMPLES_RANGE"][lignol][0],
                config["MIN_SAMPLES_RANGE"][lignol][1] + 1,
            )
        ],

# Find the best value for the HDBSCAN min_samples parameter.
rule min_samples:
    input:
        models=lambda wildcards: [
            str(
                Path(
                    analysis_dir,
                    "clustering/HDBSCAN_model",
                    "nclusters_" + str(min_samples) + ".dat",
                )
            )
            for min_samples in range(
                config["MIN_SAMPLES_RANGE"][wildcards.lignol][0],
                config["MIN_SAMPLES_RANGE"][wildcards.lignol][1] + 1,
            )
        ],
        script=Path("../scripts/min_samples.py"),
    output:
        html=str(analysis_dir / "clustering/min_samples.html"),
        png=str(analysis_dir / "clustering/min_samples.png"),
        png_writing = str(writing_dir / "fig-nclusters_{lignol}" / "min_samples.png"),
        min_samples=str(analysis_dir / "clustering/min_samples_best.dat"),
    params:
        drctry=str(analysis_dir / "clustering"),
        nclusters=lambda wildcards: config["NCLUSTERS_TARGET"][wildcards.lignol],
        min_samples_min=lambda wildcards: config["MIN_SAMPLES_RANGE"][
            wildcards.lignol
        ][0],
        min_samples_max=lambda wildcards: config["MIN_SAMPLES_RANGE"][
            wildcards.lignol
        ][1],
    shell:
        """
        python {input.script} --lignol {wildcards.lignol} \
            --min_samples_rng {params.min_samples_min} {params.min_samples_max} \
            --dir {params.drctry}
        cp {output.png} {output.png_writing}
        """

rule min_sampless:
    input:
        expand(
            rules.min_samples.output,
            lignol=config["LIGNOLS"],
        ),

# Save the best HDBSCAN model.
rule best_model:
    input:
        colvar=rules.plumed_CVs.output.colvar,
        script=Path("../scripts/clustering.py"),
        read_script=Path("../scripts/read_colvar_file.py"),
        min_samples=rules.min_samples.output.min_samples,
    output:
        model=str(analysis_dir / "clustering/HDBSCAN_model/best_model.joblib"),
    params:
        min_samples=lambda wildcards: np.loadtxt(f"../analysis/{wildcards.lignol}/one_BCD/clustering/min_samples_best.dat", dtype=int)
    shell:
        "python {input.script} --lignol {wildcards.lignol} --colvar {input.colvar} --min_samples {params.min_samples} --model_file {output.model}"

rule best_models:
    input:
        expand(rules.best_model.output, lignol=config["LIGNOLS"]),   

# Plot the clusters.
rule plot_clusters:
    input:
        script=Path("../scripts/plot_clusters.py"),
        model=rules.best_model.output.model,
        colvar=rules.plumed_CVs.output.colvar,
    output:
        html=str(analysis_dir / "clustering/clusters.html"),
        json=str(analysis_dir / "clustering/clusters.json"),
        json_writing=str(writing_dir / "fig-clusters_{lignol}" / "clusters.json"),
        png=str(analysis_dir / "clustering/clusters.png"),
        cluster_labels=str(analysis_dir / "clustering/cluster_labels.dat"),
        nclusters=str(analysis_dir / "clustering/nclusters.dat"),
    params:
        drctry=str(analysis_dir / "clustering"),
    shell:
        """
        python ../scripts/plot_clusters.py \
            --lignol {wildcards.lignol} \
            --colvar {input.colvar} \
            --model {input.model} \
            --dir {params.drctry}
        cp {output.json} {output.json_writing}
        """

rule plot_clusterss:
    input:
        expand(
            rules.plot_clusters.output, lignol=config["LIGNOLS"]
        ),

# Save the configurations belonging to each cluster to separate XTC files.
rule cluster_configs:
    input:
        cluster_labels=rules.plot_clusters.output.cluster_labels,
        script=Path("../scripts/cluster_configs.py"),
        gro=str(output_dir / "nowater.gro"),
        xtc=str(output_dir / "nowater.xtc"),
        nclusters=rules.plot_clusters.output.nclusters,
    output:
        xtc=str(analysis_dir / "clustering/cluster_configs_{cluster_id}.xtc"),
    shell:
        "python {input.script} --labels {input.cluster_labels} --gro {input.gro} --xtc_in {input.xtc} --xtc_out {output.xtc} --cluster_id {wildcards.cluster_id}"


rule cluster_configss:
    input:
        lambda wildcards: [
            str(
                Path(
                    "../analysis",
                    lignol,
                    "one_BCD/clustering",
                    "cluster_configs_" + str(cluster_id) + ".xtc",
                )
            )
            for lignol in config["LIGNOLS"]
            for cluster_id in range(
                int(
                    np.loadtxt(
                        str(
                            Path(
                                "../analysis",
                                lignol,
                                "one_BCD/clustering/nclusters.dat",
                            )
                        )
                    )
                )
            )
        ],


# Draw arrows from lignin head to tail, and from the cyclodextrin COM to the cyclodextrin secondary face oxygen atoms onto a snapshot in VMD.
rule vmd_arrows:
    input:
        plumed=str(rules.plumed_CVs.input.plumed),
        script=Path("../scripts/vmd_arrows.py"),
    output:
        tcl=str(analysis_dir / "arrows.tcl"),
    params:
        arrow_length=config["VMD_ARROW_LENGTH"],
    shell:
        "python {input.script} --plumed {input.plumed} --tcl {output.tcl} --arrow_length {params.arrow_length}"


rule vmd_arrowss:
    input:
        expand(rules.vmd_arrows.output, lignol=config["LIGNOLS"]),


# Add example configurations with vector arrows to the cluster plots.
rule cluster_configs_to_cluster_plots:
    input:
        cluster_plot=rules.plot_clusters.output.png,
        cluster_configs=[
            str(file)
            for file in Path(
                "../analysis/{lignol}/one_BCD/clustering"
            ).glob("cluster_configs_*.tga")
        ],
        script=Path("../scripts/process_cluster_config_images.sh"),
    output:
        cluster_configs_plot=Path(
            "../analysis/{lignol}/one_BCD/clustering/clusters_configs.png"
        ),
        cluster_configs_plot_writing=str(writing_dir / "fig-clusters_{lignol}" / "clusters_configs.png"),
    shell:
        """
        bash {input.script} {wildcards.lignol}
        cp {output.cluster_configs_plot} {output.cluster_configs_plot_writing}
        """


rule cluster_configs_to_cluster_plotss:
    input:
        expand(
            rules.cluster_configs_to_cluster_plots.output,
            lignol=config["LIGNOLS"],
        ),
