
from os import listdir
from os.path import isfile, join
import xml.etree.cElementTree as ET
from myKeggPathway import myKeggPathway
import random

rsample = random.Random();

pathway_path = '../data/pathway with action'
pathwayfiles = [ f for f in listdir(pathway_path) if isfile(join(pathway_path,f)) ]

#store the probability of all genesC:\Users\Mansoori\Desktop\fatima proposal\implementation\parser\data\pathway with action
prob = {}
pathways = []

for f in pathwayfiles:
    xmlfile = join(pathway_path,f)
    tree = ET.parse(xmlfile)

    pathway = myKeggPathway();

    # Get pathway title (store it in pathway.title)
    pathway.title = tree.getroot().get('title')
    pathway.name = tree.getroot().get('name')
    pathway.id = tree.getroot().get('number')

    print pathway.title
    print pathway.id

    if(pathway.id =="05120" ):
        pathway.id = pathway.id;

    # parse and add nodes
    for entry in tree.getiterator('entry'):
        node_type = entry.get('type') # can be ('gene', 'compound', 'map'..)
        #get the name of the node the name is like hsa:x, for node with gene type, this name is the name of the protein
        name = entry.get('name')
        #the id of the node in pathway, this id is used for referencing the node
        node_id = entry.get('id')

        if(node_type == 'gene'):
            pathway.add_protein(protein_data={'id':node_id, 'name':name})

        if(node_type == 'compound'):
            pathway.add_compound(comp_data={'id':node_id,'name':name})

        if(node_type == 'group'):
            group = []
            for c in entry.getiterator('component'):
                group.append(c.get('id'));
            pathway.add_groups(group_data={'id':node_id, 'components':group})


    for rel in tree.getiterator('relation'):
        e1 = rel.get('entry1')
        e2 = rel.get('entry2')
        type=rel.get('type')

        t = rel.find('subtype');
        if(t != None):
            rel_type = t.get('name');
        else:
            rel_type='None'

        #edit groups
        groups = pathway.get_groups();
        if(e1 in groups):
            cmp = groups[e1]['components'];
            for c in cmp:
                pathway.add_relation(data={'type': type, 'type_name': rel_type, 'first': c, 'second': e2, 'prob': 1})
        elif(e2 in groups):
            cmp = groups[e2]['components'];
            for c in cmp:
                pathway.add_relation(data={'type': type, 'type_name': rel_type, 'first': e1, 'second': c, 'prob': 1})
        else:
            pathway.add_relation(data={'type': type, 'type_name': rel_type, 'first': e1, 'second': e2, 'prob': 1})

    for act in tree.getiterator('action'):
        type = act.get('type');
        ent1 = act.get('entry1');
        pathway.add_action(data={'type':type, 'entry1':ent1})

    pathways.append(pathway);


#for each relation in pathways calculate the probabilities
for pathway in pathways:
    rels = pathway.get_relations();

    #remove componet relations
    relations = pathway.get_relations()
    components = pathway.get_compounds();
    new_relations = []
    for r in relations:
        rdic = relations[r];
        for idx in rdic:
            rel = rdic[idx];
            if (rel['type'] == 'PCrel'):
                first = rel['first'];
                if (first in components):
                    psecond = rel['second']
                    if (first in relations):
                        rcmp = relations[first];
                        for rc in rcmp:
                            pfirst = rcmp[rc]['first'];
                            new_relations.append({'type': 'PPrel', 'type_name': rel['type_name'], 'first': pfirst, 'second': psecond, 'prob': 1})

    for nr in new_relations:
        pathway.add_relation(nr);

    poplist = [];
    for r in relations:
        if r in components:
            poplist.append(r);
    for pl in poplist:
        pathway.remove_relation(pl);


    poplist=[]
    for r in relations:
        rdic = relations[r];
        temp = [idx for idx in rdic if rdic[idx]['type'] != 'PCrel']
        if(len(temp)!=0):
            rtemp={}
            for i in xrange(len(temp)):
                rtemp[i] = rdic[temp[i]];
            pathway.add_relaions_of_protein(r, rtemp)
        else:
            poplist.append(r);

    for pl in poplist:
        pathway.remove_relation(pl);



    pcnt = 0;
    rcnt = 0;
    relArrKeys = rels.keys()
    for relKey in relArrKeys:
        e2 = relKey;
        temprel = []
        relArr = rels[e2];
        for t in relArr:
           temprel.append(t);
        for rel in temprel:
            if(relArr[rel]['type_name']!='identiacal'):
                e1 = relArr[rel]['first'];
                #get the name of the protein e1
                if (e1 in pathway.get_proteins() and e2 in pathway.get_proteins()):
                    e1name = pathway.get_proteins()[e1]['name'];
                    e2name = pathway.get_proteins()[e2]['name'];
                    #search other pathways for these genes
                    for path in pathways:
                        exist = False;
                        proteins = path.get_proteins();
                        for p1 in proteins:
                            if(proteins[p1]['name'] == e1name):
                                for p2 in proteins :
                                    if(proteins[p2]['name']==e2name):
                                        pcnt+=1.0;
                                        pp1 = p1;
                                        pp2 = p2;
                                        exist = True;
                                if(exist):
                                    #search the relation to see if they have relation to each other
                                    relations = path.get_relations();
                                    if(relations.has_key(pp2)):
                                        for i in relations[pp2]:
                                            if relations[pp2][i]['first']==pp1:
                                                rcnt+=1.0;

                    if(e1=='14' and e2=='4'):
                        e1 = e1;

                    if (rcnt==0):
                        rcnt=1.0;
                    if(pcnt==0):
                        pcnt=1.0;

                    prob = rcnt/pcnt;
                    if(prob==1.0):
                        prob = 0.95

                    pathway.add_relation({'type':relArr[rel]['type'], 'type_name':relArr[rel]['type_name'], 'first':e1, 'second':e2, 'prob': prob})
                    rcnt = 0;
                    pcnt = 0;

                else:
                    pathway.add_relation(
                        {'type': relArr[rel]['type'], 'type_name': relArr[rel]['type_name'], 'first': e1, 'second': e2,
                         'prob': 0.95})
            else:
                #identical relation
                pathway.add_relation({'type':relArr[rel]['type'], 'type_name':relArr[rel]['type_name'], 'first':e1, 'second':e2, 'prob': 0.95})

#write pathway to files.
#write proteins.
for pathway in pathways:

    fout = open('../data/pathway with action/with prob/hsa'+pathway.id+'.xml', 'w')

    fout.write('<?xml version="1.0"?>\n')
    fout.write('<pathway name="'+pathway.name+'" org="hsa" '+' number="'+pathway.id+'" title="'+pathway.title+'">\n');

    proteins = pathway.get_proteins();
    for p in proteins:
        fout.write('<entry id="'+proteins[p]['id']+'" type="protein" name="'+proteins[p]['name']+'"/>\n')

    compounds = pathway.get_compounds();
    for c in compounds:
        fout.write('<entry id="'+compounds[c]['id']+'" type="compound" name="'+compounds[c]['name']+'"/>\n');

    groups = pathway.get_groups();
    for g in groups:
        fout.write('<entry id="'+groups[g]['id']+'" type="group">\n');
        comp = groups[g]['components'];
        for c in comp:
            fout.write('<component id="'+c+'"/>\n');
        fout.write('</entry>\n')


    relations = pathway.get_relations()
    components = pathway.get_compounds();
    for r in relations:
        rdic = relations[r];
        for idx in rdic:
            rel = rdic[idx];
            fout.write(
                    '<relation entry1="' + rel['first'] + '" entry2="' + rel['second'] + '" type="' + rel['type'] +
                    '" prob="' + str(round(rel['prob'], 10)) + '">\n');

            fout.write('<subtype name="' + rel['type_name'] + '"/>\n');
            fout.write('</relation>\n');

    #write identical relations
    relations = pathway.get_relations()
    for r in relations:
        rdic = relations[r];
        for idx in rdic:
            rel = rdic[idx];
            if(rel['type_name']=='identical'):
                fout.write('<relation entry1="'+rel['first']+'" entry2="'+rel['second']+'" type="'+rel['type']+
                        '" prob="'+str(round(rel['prob'],10))+'">\n');

                fout.write('<subtype name="'+rel['type_name']+'"/>\n');
                fout.write('</relation>\n');

    #write not identical relations
    relations = pathway.get_relations()
    for r in relations:
        rdic = relations[r];
        for idx in rdic:
            rel = rdic[idx];
            if(rel['type_name']!= 'identical'):
                #fout.write('<relation entry1="'+rel['first']+'" entry2="'+rel['second']+'" type="'+rel['type']+
                #        '" prob="'+str(round(rel['prob'],10))+'">\n');
                fout.write(
                    '<relation entry1="' + rel['first'] + '" entry2="' + rel['second'] + '" type="' + rel['type'] +
                          '" prob="'+str(round(rsample.uniform(0,1),10))+'">\n');

                fout.write('<subtype name="'+rel['type_name']+'"/>\n');
                fout.write('</relation>\n');

    #write actions
    actions = pathway.get_actions();
    for acttype in actions:
        acts = actions[acttype];
        for act in acts:
            fout.write('<action type="'+acts[act]['type']+'" entry1="'+acts[act]['entry1']+'" />\n');

    fout.write('</pathway>');






