import os

configfile: 'plumed.json'

# Need to get masses into pdb occupancy column for PLUMED.
# Note that masses guessed by MDAnalysis do not match masses from topology exactly.
rule nowater_pdb:
    input:
        gro = os.path.join('..', 'output', '{concentration}M_NaCl', '{lignol}', 'one_BCD', 
            'nowater.gro')
    output:
        pdb = os.path.join('..', 'analysis', '{concentration}M_NaCl', '{lignol}', 'nowater.pdb')
    run:
        import MDAnalysis as mda
        u = mda.Universe(input.gro, input.gro)
        u.add_TopologyAttr('occupancies')
        u.atoms.occupancies = u.atoms.masses
        u.atoms.write(output.pdb)

rule nowater_pdbs:
    input:
        expand(rules.nowater_pdb.output.pdb, concentration=['0.0'], lignol=config['LIGNOLS'])

rule plumed_CVs:
    input:
        plumed = os.path.join('{concentration}M_NaCl', '{lignol}', 'plumed.inp'),
        traj = os.path.join('..', 'output', '{concentration}M_NaCl', '{lignol}', 'one_BCD', 'nowater.xtc'),
        pdb = rules.nowater_pdb.output.pdb
    output:
        colvar = os.path.join('..', 'analysis', '{concentration}M_NaCl', '{lignol}', 'one_BCD', 
            'colvar.dat')
    shell:
        'plumed driver --plumed {input.plumed} --timestep 4 --mf_xtc {input.traj} --pdb {input.pdb}'

rule plumed_CVss:
    input:
        expand(rules.plumed_CVs.output, concentration=['0.0'], lignol=config['LIGNOLS'])

rule plots_2D:
    input:
        colvar = rules.plumed_CVs.output.colvar,
        script = '../scripts/plots_2D.py'
    output:
        plot_dnorm_dtang = os.path.join('..', 'analysis', '{concentration}M_NaCl', '{lignol}', 'one_BCD', 'dnorm_dtang.png'),
        plot_dnorm_dang = os.path.join('..', 'analysis', '{concentration}M_NaCl', '{lignol}', 'one_BCD', 'dnorm_orient.png'),
        plot_dtang_orient = os.path.join('..', 'analysis', '{concentration}M_NaCl', '{lignol}', 'one_BCD', 'dtang_orient.png')
    params:
        output_dir = os.path.join('..', 'analysis', '{concentration}M_NaCl', '{lignol}', 'one_BCD')
    shell:
        'python ../scripts/plots_2D.py --lignol {wildcards.lignol} --colvar {input.colvar} --outdir {params.output_dir}'

rule plots_2Ds:
    input:
        expand(rules.plots_2D.output, concentration=['0.0'], lignol=config['LIGNOLS'])