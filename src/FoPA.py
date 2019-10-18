__author__ = 'Mansoori'



import xml.etree.cElementTree as ET
import logging
import sys
from prism_trim_parser import prism_trim_parser
from sort import sorting
from KeggPathway import KeggPathway

from os import listdir
from os.path import isfile, join

import ConfigParser


####################


####################
confige_file = './config.cfg'
###################################


#read the config file
config = ConfigParser.ConfigParser()
config.read(confige_file)
#read the config data
##################################
#pathes
FoPA_path = config.get('paths','FoPA_path')
gene_path = config.get('paths','gene_path')
pathway_path = config.get('paths','pathway_path')
null_path = config.get('paths','null_path')

freq_file = config.get('files','freq_file')
result_file = config.get('files','result_file')
ignore_file = config.get('files','ignore_file')

random_build = config.get('Random','Build')
random_number = config.get('Random','Number')


#read the arguments
#if(len())
if(len(sys.argv) > 4):
    new_continue = sys.argv[1].lstrip();
    deg_file = sys.argv[2].lstrip();
    all_file = sys.argv[3].lstrip();
    out_file = sys.argv[4].lstrip();
else:
    print("Command format: FoPA -(n/c) 'DEGs file' 'all genes file' 'result file'.")
    print("Read the ReadMe file for help.")
    sys.exit(-1)

#read freq file and build the weights
freq_raw = [line.rstrip('\n') for line in open(freq_file)]
ws={}
wmax = 0;
wmin = 20000;
for d in freq_raw:
    splitd = d.split('\t\t');
    w = abs(float(splitd[2]));
    ws[splitd[0]] = w;
    if(w < wmin):
        wmin = w;
    if(w > wmax) :
        wmax = w;
#build the weights
for we in ws:
    w = ws[we]
    ws[we] = ((wmax-w)/(wmax-wmin))**(0.5);


#read result file to continue from previous if the option is -c
donePathwaysFiles=[];
if(new_continue == 'c'):
    res_file = open(result_file,'r');
    donePathwaysRows = [line.rstrip('\n') for line in res_file];
    res_file.close();
    for donePathwaysRow in donePathwaysRows:
        splitd = donePathwaysRow.split('\t');
        donePathwaysFiles.append("hsa"+splitd[0]+".xml");


#read ig file
ignorePathwaysFiles = [];
ignorePathways = [line.rstrip('\n') for line in open(ignore_file,'r')];
for ignorePathwaysRow in ignorePathways:
    ignorePathwaysFiles.append(ignorePathwaysRow+".xml");


#Read pathways name
pathwayfiles = [ f for f in listdir(pathway_path) if isfile(join(pathway_path,f)) ]


#parse the pathway, build the model and calculate the score
for f in pathwayfiles:
    # read DEG file
    degs_raw = [line.rstrip('\n') for line in open(join(gene_path, deg_file))]
    # read all gene file
    all_raw = [line.rstrip('\n') for line in open(join(gene_path, all_file))]
    degs = {}
    for d in degs_raw:
        splitd = d.split(' ');
        degval = abs(float(splitd[2]));
        degs["hsa" + splitd[1]] = degval;
    # enter all genes t value
    all_genes = {}
    tmax = 0.0;
    tmin = 2000.0;
    for g in all_raw:
        splitg = g.split(' ');
        tval = abs(float(splitg[2]))
        all_genes["hsa" + splitg[1]] = tval;
        if (tmax < tval):
            tmax = tval;
        if (tmin > tval):
            tmin = tval;
    # normalize t value
    # we want that omit the zero so we use tmin/10000 instead of tmin
    tmin = tmin/10000;
    for g in all_genes:
        tval = all_genes[g];
        all_genes[g] = float(tval - tmin) / (tmax - tmin);


    #ignore pathways in result file or ignore file
    if ((f in donePathwaysFiles) or (f in ignorePathwaysFiles)): continue;

    #parse the pathway file
    xmlfile = join(pathway_path, f);
    tree = ET.parse(xmlfile)

    pathway = KeggPathway();
    pathway.__init__();

    # Get pathway title name and id and store
    pathway.title = tree.getroot().get('title')
    pathway.name = tree.getroot().get('name')
    pathway.id = tree.getroot().get('number')

    print(pathway.title + '\n');

    # parse pathway and add nodes
    maxnodeid = 0;
    for entry in tree.getiterator('entry'):
        node_type = entry.get('type')  # can be ('gene', 'compound', 'map'..)
        # get the name of the node the name is like hsa:x, for node with gene type, this name is the name of the protein
        name = entry.get('name')
        # for using the name for naming the variable in prism, the : is removed from it
        if (name != None):
            name = name.replace(":", "")
        # the id of the node in pathway, this id is used for referencing the node
        node_id = entry.get('id')

        if (int(node_id) > maxnodeid):
            maxnodeid = int(node_id);

        if (node_type == 'protein'):
            genes = pathway.split_genes(name)

            if (len(genes) == 1):
                pname = name;

                if(name in all_genes):
                    prob = all_genes[name]*ws[name];
                else :
                    prob = 0;
                # add protein in pathway's proteins
                pathway.add_protein(protein_data={'id': node_id, 'name': pname, 'prob': prob})
                pathway.add_protein_genes(data={'pname': pname, 'genes': genes});
            else:
                pname = genes[0]

                prob = 0;
                cntg = 0;
                for g in genes:
                    # add gene in pathway's genes
                    if(g in all_genes):
                        prob += all_genes[g]*ws[g];
                    cntg+=1;
                    pathway.add_gene(gene_data={'name': g, 'protein': pname})
                pathway.add_protein_genes(data={'pname': pname, 'genes': genes});

                prob = prob / cntg;
                # add protein in pathway's proteins
                pathway.add_protein(protein_data={'id': node_id, 'name': pname, 'prob': prob})
        elif (node_type == 'group'):
            group = []
            group_names = []
            for c in entry.getiterator('component'):
                pid = c.get('id');
                pid = pathway.get_uniq_id(pid);
                proteins = pathway.get_proteins();
                if (pid in proteins):
                    pname = pathway.get_proteins()[pid]['name'];
                else:
                    pname = pid
                if (pname not in group_names):
                    group_names.append(pname);
                    group.append(pid);
            pathway.add_groups(group_data={'id': node_id, 'components': group})

        elif (node_type == 'compound'):
            comps = pathway.split_genes(name);
            if (len(comps) == 1):
                pathway.add_compound(comp_data={'id': node_id, 'name': name});
            else:
                pathway.add_compound(comp_data={'id': node_id, 'name': comps[0]});
    for rel in tree.getiterator('relation'):
        e1 = rel.get('entry1')
        e2 = rel.get('entry2')
        type = rel.get('type')
        prob = rel.get('prob');

        if (rel.find('subtype') != None):
            rel_type = rel.find('subtype').get('name')
        else:
            rel_type = "NONE";
            # pathway.add_edge(data={'type':rel_type, 'first':e1, 'second':e2})
        if (rel_type == 'identical'):
            pathway.add_identical_relation(rel_data={'first': e1, 'second': e2});
        else:
            pathway.add_relation(data={'type': type, 'type_name': rel_type, 'first': e1, 'second': e2, 'prob': prob})

    #compute the probability
    myprism_trim_parser = prism_trim_parser();
    myprism_trim_parser.__init__();
    path_prob = myprism_trim_parser.calculate_score(pathway, degs, all_genes);
    print '\n',
    print path_prob;
    outputf = open(result_file, 'a');
    outputf.write(pathway.id + '\t' + pathway.title + '\t' + str(path_prob) + '\n');
    outputf.close();


    # iterate some time for each pathway and permute the labels
    proteins = pathway.get_proteins();
    degnames = degs.keys();
    degcnt = 0;
    protein_genes = pathway.get_protein_genes();
    degs_in_pathway = [];
    all_pathway_genes = [];
    for p in proteins:
        pname = proteins[p]['name'];

        genes = protein_genes[pname];

        isDEG = False;
        for gene in genes:
            if (gene not in all_pathway_genes):
                all_pathway_genes.append(gene);
                if gene in degnames:
                    degs_in_pathway.append(gene);
                    isDEG = True;

    if(random_build == 'ALL'):
        for s in xrange(random_number):
            print s;
            # read DEG file
            degs_raw = [line.rstrip('\n') for line in open(join(gene_path, deg_file))]
            # read all gene file
            all_raw = [line.rstrip('\n') for line in open(join(gene_path, all_file))]

            degs = {}
            for d in degs_raw:
                splitd = d.split(' ');
                degval = abs(float(splitd[2]));
                degs["hsa" + splitd[1]] = degval;

            # enter all genes t value
            all_genes = {}
            tmax = 0.0;
            tmin = 2000.0;
            for g in all_raw:
                splitg = g.split(' ');
                tval = abs(float(splitg[2]))
                all_genes["hsa" + splitg[1]] = tval;
                if (tmax < tval):
                    tmax = tval;
                if (tmin > tval):
                    tmin = tval;
            # normalize t value
            tmin = tmin / 10000;
            for g in all_genes:
                tval = all_genes[g];
                all_genes[g] = float(tval - tmin) / (tmax - tmin);

            #update proteins prob
            for p in proteins:
                pname = proteins[p]['name'];
                genes = protein_genes[pname];
                prob = 0;
                gcnt = 0;

                if (len(genes) == 1):
                    if (pname in all_genes):
                        prob = all_genes[pname] * ws[pname];
                else:

                    prob = 0;
                    cntg = 0;
                    for g in genes:
                        if(g in all_genes):
                            prob += all_genes[g] * ws[g];
                        cntg += 1;
                    prob = prob / cntg;


                pathway.update_protein_prob(p, prob);


            myprism_trim_parser = prism_trim_parser();
            myprism_trim_parser.__init__();
            path_prob = myprism_trim_parser.calculate_prob(pathway, degs, all_genes);
            print '\n'
            print path_prob;
            randomf = open('./data/temp/p'+pathway.id+'.txt', 'w');
            randomf.write(pathway.id + '\t' + pathway.title + '\t' + str(path_prob) + '\n');
            randomf.close();


if (random_build == 'ALL'):
            mysorting = sorting();
            mysorting.sort(out_file)
else:
    mysorting = sorting();
    mysorting.sort_stored_random(out_file)
