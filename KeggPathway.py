__author__ = 'Mansoori'

import string

class KeggPathway():

    title = '';
    name='';
    id='';

    #store all nodes of pathway (genes, proteins, active proteins)
    nodes = {};
    #store genes of the pathway
    genes = {};
    #store proteins of the pathway
    proteins = {};
    proteins_names = {};
    proteins_id_map = {};

    protein_genes = {}

    compounds = {};


    #store relations of the pathway
    relations = {};

    reverse_relations = {};

    groups = {};

    actions = {};

    identical_relation = {};

    def __init__(self):

        self.title = '';
        self.name = '';
        self.id='';

        self.proteins = {};
        self.proteins_names = {};
        self.proteins_id_map = {};

        #store all nodes of pathway (genes, proteins, active proteins)
        self.nodes = {};
        #store genes of the pathway
        self.genes = {};


        self.protein_genes = {}
        self.compounds = {};

        self.groups = {};

        #store relations of the pathway
        self.relations = {};

        self.reverse_relations = {};

        self.actions ={};

    def add_node(self, node_data=None):
        self.nodes[node_data['name']]=node_data;

    def add_gene(self, gene_data=None):
        self.genes[gene_data['name']] = gene_data;

    def add_protein_genes(self, data=None):
        self.protein_genes[data['pname']] = data['genes'];

    def add_protein(self, protein_data=None):
        self.proteins[protein_data['id']]= protein_data;
        self.add_uniq_protein(protein_data);

    def update_protein_prob(self, proteind_id, prob):
        self.proteins[proteind_id]['prob'] = prob;

    def add_compound(self, comp_data=None):
        self.compounds[comp_data['id']]=comp_data;

    def add_groups(self, group_data=None):
        self.groups[group_data['id']] = group_data;

    def add_identical_relation(self, rel_data=None):
        self.identical_relation[rel_data['first']] = rel_data;

    def remove_relation(self, second):
        self.relations.pop(second);


    def add_relaions_of_protein(self, second, rels):
        self.relations[second] = rels;


    def add_relation(self, data=None):

        if(self.proteins.has_key(data['first'])):
            data['first'] = self.proteins_id_map[data['first']];
        if(self.proteins_id_map.has_key(data['second'])):
            data['second'] = self.proteins_id_map[data['second']];

        #check for identical relation
        if(self.identical_relation.has_key(data['first'])):
            data['first'] = self.identical_relation[data['first']]['second'];
        if(self.identical_relation.has_key(data['second'])):
            data['second'] = self.identical_relation[data['second']]['second'];

        if(self.relations.has_key(data['second'])):
            exist = False;
            for i in self.relations[data['second']] :
                if self.relations[data['second']][i]['first']==data['first']:
                    self.relations[data['second']][i] = data;
                    exist = True
            if exist==False:
                d = len(self.relations[data['second']])
                self.relations[data['second']][d] = data;
        else:
            dict = {}
            dict[0] = data
            self.relations[data['second']] = dict;


        #also add reverse relation
        if(self.reverse_relations.has_key(data['first'])):
            exist = False;
            for i in self.reverse_relations[data['first']] :
                if self.reverse_relations[data['first']][i]['second']==data['second']:
                    self.reverse_relations[data['first']][i] = data;
                    exist = True
            if exist==False:
                d = len(self.reverse_relations[data['first']])
                self.reverse_relations[data['first']][d] = data;
        else:
            dict = {}
            dict[0] = data
            self.reverse_relations[data['first']] = dict;




    def add_action(self, data=None):

        if(data['entry1'] in self.proteins_id_map):
            data['entry1'] = self.proteins_id_map[data['entry1']];


        if(self.actions.has_key(data['type'])):
            d = len(self.actions[data['type']]);
            self.actions[data['type']][d] = data
        else:
            dict = {}
            dict[0] = data;
            self.actions[data['type']] = dict;

    #remove proteins with redundant name
    def add_uniq_protein(self, protein_data = None):
        pname = protein_data['name'];
        pid = protein_data['id'];
        if (self.proteins_names.has_key(pname)):
            self.proteins_id_map[pid] = self.proteins_names[pname];
        else:
            self.proteins_names[pname] = pid;
            self.proteins_id_map[pid] = pid;


    def get_uniq_id(self, id):

        return self.proteins_id_map[id];

    def get_protein_genes(self):
        return self.protein_genes;

    def get_genes(self):
        return self.genes;

    def get_proteins(self):
        return self.proteins;

    def get_relations(self):
        return self.relations;

    def get_nodes(self):
        return self.nodes;

    def get_compounds(self):
        return self.compounds;

    def get_groups(self):
        return self.groups;

    def get_actions(self):
        return self.actions;

    def get_identical_relations(self):
        return self.identical_relation;

    def get_rrelations(self):
        return self.reverse_relations;




    '''def get_genes(self):
        genes = {}
        id = 0
        for n in self.nodes:
            if self.nodes[n]['type'] == 'gene':
                name = self.nodes[n]['name']
                splitName = name.split(' ')
                for s in splitName:
                    genes[id] = s
                    id = id + 1

        return genes
    '''

    def split_genes(self, nodename):
        genes = []
        splitName = nodename.split(' ')
        for s in splitName:
            genes.append(s)

        return genes


    def add_edge(self,data=None):
        self.edge.append(data)


    def get_node(self, n):
        return self.nodes[n]



