__author__ = 'Mansoori'


from subprocess import call
from prism_trim import prism_trim

import ConfigParser

confige_file = './config.cfg'


class prism_trim_parser():

    prism_trim;

    prism_path = '.\data\prism output';

    pathway_prob = [];
    all_gene_pathway_prob = {}

    def __init__(self):
        self.prism_trim = prism_trim();
        self.pathway_prob = [];

    def var_print(self, vars=[], svar=None, firsts=[], relations={}):

        if(svar not in vars):
            vars.append(svar);

            seconds = relations.keys();
            for fvar in firsts:
                if (fvar in seconds) and (fvar not in vars) :
                    fs = relations[fvar];
                    new_firsts = [];
                    for fsi in fs :
                        new_firsts.append(fs[fsi]['first']);
                    self.var_print(vars, fvar, new_firsts, relations);
        for fvar in firsts :
            if (fvar not in vars):
                vars.append(fvar);


    def find_path(self,root, pvars=[], relations={}):
        if(root not in pvars):
            pvars.append(root);
            seconds = relations.keys();
            if(root in seconds):
                root_relations = relations[root];
                for fri in root_relations:
                    fr = root_relations[fri]['first'];
                    if (fr not in pvars):
                        self.find_path(fr, pvars, relations);


    def remove_unwanted_vars(self, all_vars=[], firsts={}, seconds={}):

        for k in seconds.keys():
            if(k not in all_vars):
                ss = seconds[k];
                for s in ss:
                    fi = ss[s]['first'];

                    if (fi in firsts):
                        ff = firsts[fi];
                        ks = ff.keys();
                        for cnt in ks:
                            if(ff[cnt]['second']==k):
                                firsts[fi].pop(cnt);

                    if(fi in firsts and (firsts[fi]) == {}) :
                        firsts.pop(fi);

                seconds.pop(k);


    def calculate_score(self, pathway=None, degs=None, all_genes=None):

        relations = pathway.get_relations();
        proteins = pathway.get_proteins();
        components = pathway.get_compounds();
        actions = pathway.get_actions();
        groups = pathway.get_groups();

        #free pathway probability
        self.pathway_prob = [];

        #find pathway last points
        pathway_last_points = [];
        first_of_pathways = pathway.get_rrelations();
        seconds_of_pathways = pathway.get_relations();
        for s in seconds_of_pathways.keys():
            if (s not in first_of_pathways.keys()) :
                pathway_last_points.append(s);


        # add variable
        all_vars = [];
        for s in seconds_of_pathways.keys():
            all_vars.append(s);
        for f in first_of_pathways.keys():
            if f not in all_vars:
                all_vars.append(f);
        old_last_points = {}
        #copy pathway relations
        firsts = {}
        firsts_temp = pathway.get_rrelations().copy();
        for f in firsts_temp:
            l = {}
            ll = firsts_temp [f];
            for g in ll:
                l[g] = ll[g];
            firsts[f] = l;

        #copy pathway relations
        seconds = {}
        seconds_temp = pathway.get_relations().copy();
        for s in seconds_temp:
            l = {}
            ll = seconds_temp[s];
            for g in ll:
                l[g] = ll[g];
            seconds[s] = l;

        self.prism_trim.set_firsts_seconds(firsts, seconds);
        self.prism_trim.set_all_variables(all_vars);

        while (True):

            trim_points = [];
            begins = {}
            ppoints = [];
            new_last_points = {};

            self.prism_trim.trim(old_last_points, begins, trim_points, new_last_points);

            if  len(trim_points) == 0 :
                break;
            else:
                #after trim, the trim points are the variables that should be in the prism file
                #the probability of the last points should be calculated.
                self.calcute_trim_prob(old_last_points, new_last_points, trim_points, pathway, degs, all_genes);
                old_last_points.clear();
                for n in new_last_points:
                    old_last_points[n] = new_last_points[n];
                    if(n in pathway_last_points):
                        self.pathway_prob.append(new_last_points[n]['prob4']);
            print '.',

        sum_prob = 0.0;
        for pp in self.pathway_prob:
            sum_prob += float(pp);
            #if float(pp) > max :
            #    max = float(pp);

        return sum_prob;


    def run_prism_file_in_linux(self, filename=None, profile=None):

        path = '/home/ahmadreza/fatima/prism-4.3-src/bin/data';
        #ppath = 'C:/Program Files/prism-4.1.beta2/bin';

        #call([ppath+"\\prism "+path+'\\'+filename+' '+path+'\\prop.csl >out']);
        #f = open(ppath+'/test1.sm','r');

        #call(["myprism.bat",path+'/hsa03320.sm"',path+'/prop.csl"','>'+path+'/out"']);

        #true call
        call(["prism", path+"/"+filename, path+'/'+profile, '>/out' ]);


        f = open('/out', 'r');

        for line in f.readlines():
            if(line.find('Result:')!=-1):
                s = line.split(' ');
                return s[1];




    def run_prism_file_in_windows(self, filename=None, profile=None):

        # read the config file
        config = ConfigParser.ConfigParser()
        config.read(confige_file)
        FoPA_path = config.get('paths', 'FoPA_path')
        prism_result = config.get('paths', 'prism_result')
        path = FoPA_path+'/data/prism output';

        #call(["myprism.bat", path+'/'+filename, path+'/'+profile, '>'+prism_result+'/out' ]);
        call(["myprism.bat", './data/prism output/' + filename, './data/prism output/' + profile, '>'+prism_result+'/out']);
        f = open(prism_result+'/out', 'r');

        for line in f.readlines():
            if(line.find('Result:')!=-1):
                s = line.split(' ');
                return s[1];

    def build_prism_model(self, pathway=None, all_vars=[], last_points={}, degs=None, all_genes=None):

        relations = pathway.get_relations();
        proteins = pathway.get_proteins();
        components = pathway.get_compounds();
        actions = pathway.get_actions();
        groups = pathway.get_groups();

        filename = 'hsa'+pathway.id+'.sm';

        f = open(self.prism_path+'/'+filename, 'w')

        f.write('dtmc\n')
        f.write('module '+"hsa"+pathway.id+'\n');

        degnodes = degs.keys()
        all_vars_names = []
        protein_genes = pathway.get_protein_genes();
        for vid in all_vars:
            if vid in proteins:
                vname = proteins[vid]['name'];
                #get the genes of the vname
                genes = []
                if vname in protein_genes:
                    genes = protein_genes[vname];
                else:
                    genes.append(vname);

                if vname not in all_vars_names:
                    all_vars_names.append(vname);
                    isDEG = False;
                    isAll = False;
                    for gene in genes:
                        if gene in degnodes :
                            isDEG = True;
                        elif gene in all_genes.keys():
                            isAll = True;
                    if (isDEG==True) :
                        f.write(vname+': [-1..4] init -1; \n');
                    elif (isAll == True) :
                        f.write( vname+': [-1..4] init -1; \n');
                    else:
                        f.write( vname+': [-1..4] init -1; \n');
            elif vid in components :
                vname = components[vid]['name'];
                if(vname not in all_vars_names):
                    all_vars_names.append(vname);
                    f.write( components[vid]['name']+': [0..4] init 2; \n');
            elif vid in groups:
                for comp in groups[vid]['components']:
                    if comp in proteins:
                        vname = proteins[comp]['name'];
                        if(vname not in all_vars_names):
                            all_vars_names.append(vname);
                            if vname in degs.keys() :
                                f.write(vname+': [-1..4] init -1; \n');
                            elif vname in all_genes.keys() :
                                f.write( vname+': [-1..4] init -1; \n');
                            else:
                                f.write( vname+': [-1..4] init -1; \n');
            elif vid in self.need_action.values():
                #vname = self.need_action[vid];
                all_vars_names.append(vid);
                f.write(vid+':[0..4] init 0;\n');

        #for each variable add tree state
        for vid in all_vars:
            if vid in proteins:
                vname = proteins[vid]['name'];
                genes = []
                if vname in protein_genes:
                    genes = protein_genes[vname];
                else:
                    genes.append(vname);
                isDEG = False;
                isAll = False;
                for gene in genes:
                    if gene in degnodes:
                        isDEG = True;
                    elif gene in all_genes.keys() :
                        isAll = True;

                if (vid == '16'):
                    vid = vid;
                if(vid in last_points.keys()):
                    prob3 = float(last_points[vid]['prob3']);
                    prob4 = float(last_points[vid]['prob4']) ;
                    rprob = 1-prob3-prob4;
                    f.write('[] ' + vname + '=-1 ->' + str(prob3) + ':(' + vname + "'=3)+" + str(prob4) + ":(" + vname + "'=4)+" +
                            str(rprob) + ":(" + vname + "'=0);\n")

                else:
                    if (isDEG == True):
                        f.write('[] ' + vname + '=-1 -> 0.05:(' + vname + "'=1)+0.95:(" + vname + "'=2);\n")
                    elif isAll == True :
                        f.write('[] '+vname+'=-1 -> 0.05:('+vname+"'=0)+0.05:("+vname+"'=2)+0.9:("+vname+"'=1);\n")
                    else :
                        f.write(
                            '[] ' + vname + '=-1 -> 0.05:(' + vname + "'=1)+0.95:(" + vname + "'=0);\n")


        #find the inhibition relations
        inhibits = {};
        first = True;
        for r in relations:
            c = relations[r];
            inhibit = ""
            for k in c:
                el = c[k]
                if(el['type_name'] == 'inhibition'):
                    id = el['first'];

                    if(id in proteins):
                        p = proteins[id];
                    elif(id in components):
                        p = components[id];
                        #TO DO: search in groups
                    else:
                        p = None;

                    if(p != None):
                        if(first):
                            inhibit = "(" + p['name'] + '<3 &' + p['name'] + '> -1';
                            first = False
                        else :
                            inhibit = inhibit + '&' + p['name']+'<3 &' + p['name'] + ' > -1';
            if(first != True):
                inhibit = inhibit + ')';
                first = True
                #r indicate the second part of a relation
            inhibits[r] = inhibit;


        winhibit = {}
        firstsp = {}
        secondsp = {}
        #for each activation add a prism command
        for r in relations:
            c = relations[r];
            for k in c:
                el = c[k]
                if(el['type']=='PPrel' or el['type']=='GErel'):
                    if(el['type_name'] == 'activation' or el['type_name']=='expression' or el['type_name'] == 'binding/association'
                       or el['type_name']=='indirect effect' or el['type_name'] == 'missing interaction' or el['type_name']=='compound' or
                               el['type_name'] == 'phosphorylation'):

                        idf = el['first'];
                        ids = el['second'];

                        if((idf in proteins) and (ids in proteins) and (idf in all_vars) and (ids in all_vars)):
                            pf = proteins[idf];
                            ps = proteins[ids];

                            firstsp[idf] = pf;
                            secondsp[ids] = ps;

                            #calulate the probability of the relation. the protein probability multiply the relation probability
                            prob = float(el['prob']) * float(pf['prob']);
                            probstr = '(('+pf['name']+'+'+ps['name']+'-3)/6)*'+str(prob);
                            st = '[] '+ '('+pf['name']+'=3 | '+pf['name']+'=4)'+' & '+ '('+ps['name']+'=1 | '+ps['name']+'=2)';

                            if(inhibits[r] != ''):
                                winhibit[r] = True;
                                st = st + ' & '+inhibits[r];
                            f.write(st + ' -> ' + probstr + ':(' + ps['name'] + '\'=((' + pf['name'] + '=3&' + ps['name'] + '=1)?3:4)) +' + '1-(' + probstr + ')' + ':(' + ps['name'] + '\'= 0' + ');\n');

                        elif(((idf in groups) or (ids in groups))and(idf in all_vars) and (ids in all_vars)):
                            stn1 = '';
                            stn2 = '' ;
                            stn3 = '';
                            stsn = '';
                            stfn = '';
                            probf = 0;
                            probs = 0;
                            if(idf in groups and len(groups[idf]['components'])>1):
                                firstsp[idf] = groups[idf];
                                first = True;
                                stf = '(';
                                stfn = 'max('
                                for comp in groups[idf]['components']:
                                    if comp in proteins :
                                        pgf = proteins[comp];
                                        if (first==True) :
                                            stf += pgf['name'] + '>2 ';
                                            stfn += pgf['name'];
                                            first = False;
                                        else:
                                            stf += '& ' + pgf['name'] + '>2';
                                            stfn += ',' + pgf['name'];
                                        probf = probf+float(pgf['prob']);

                                stf += ')';
                                stfn += ')=3';
                                probf = probf / 6*len(groups[idf]['components']);

                            elif(idf in proteins):
                                pf = proteins[idf];
                                firstsp[idf] = pf;
                                probf = float(pf['prob']);
                                stf = '('+pf['name'] + '>2)';
                                stfn = pf['name'] + '=3';
                            elif(len(groups[idf]['components'])==1):
                                pid  = groups[idf]['components'][0];
                                pf = proteins[pid];
                                firstsp[idf] = pf;
                                probf = float(pf['prob']);
                                stf = '('+pf['name'] + '>2)';
                                stfn = pf['name'] + '=3';

                            sts = '';
                            if(ids in groups and len(groups[ids]['components'])>1):
                                secondsp[ids] = groups[ids];
                                first = True;
                                sts = '(';
                                stsn = 'max(';
                                for comp in groups[ids]['components']:
                                    if comp in proteins:
                                        if(first==True):
                                            first = False;
                                        else:
                                            sts += '& ';
                                            stsn += ','
                                        pgs = proteins[comp];
                                        sts += '('+ pgs['name'] + '=1 |'+pgs['name']+'=2)';
                                        stsn += pgs['name'];
                                        probs += pgs['prob'];
                                sts += ')';
                                stsn +=')=1';

                                first = True;
                                stn1 = ''
                                stn2 = ''
                                stn3 = ''
                                for comp in groups[ids]['components']:
                                    if comp in proteins:
                                        if(first == True):
                                            first = False;
                                        else:
                                            stn1 += '&' ;
                                            stn2 += '&';
                                            stn3 += '&';
                                        stn1 += '('+proteins[comp]['name']+'\'=(('+stfn+'&'+stsn+')?3:4)) ';
                                        stn2 += '('+proteins[comp]['name']+'\'= 0'+')';
                                        stn3 += '(' + proteins[comp]['name'] + '\'=0' + ')';
                            elif(ids in proteins):
                                ps = proteins[ids];
                                secondsp[ids] = ps;
                                probs = probs+float(ps['prob']);
                                sts = '('+ps['name'] +'=1 |'+ps['name']+'=2)';
                                stsn = ps['name'] +'=1';
                                stn1 = '('+ps['name'] + '\'=(('+stfn+'&'+stsn+')?3:4))';
                                stn2 = '('+ps['name'] + '\'=0)';
                                stn3 = '(' + ps['name'] + '\'=0' + ')';
                            elif(ids in groups):
                                if(len(groups[ids]['components'])==1):
                                    pid = groups[ids]['components'][0];
                                    ps = proteins[pid];
                                    secondsp[ids] = ps;

                                    sts = '('+ps['name'] +'=1 |'+ps['name']+'=2)';
                                    stsn = ps['name'] +'=1';
                                    stn1 = '('+ps['name'] + '\'=(('+stfn+'&'+stsn+')?3:4))';
                                    stn3 = '(' + ps['name'] + '\'=0' +')';
                            prob = probf * float(el['prob']);
                            f.write('[]' + stf + '&' + sts + ' -> ' + str(prob) + ':' + stn1 + ' + 1-('+str(prob)+')'+':'+stn3+';\n');

                    elif(el['type_name'] == 'inhibition' and len(c)==1):
                        idf = el['first'];
                        ids = el['second'];
                        if((idf in proteins) and (ids in proteins) and (idf in all_vars) and (ids in all_vars)):
                            pf = proteins[idf];
                            ps = proteins[ids];
                            firstsp[idf] = pf;
                            secondsp[ids] = ps;
                            prob = float(el['prob']) * float(pf['prob']);
                            probstr = '((' + pf['name'] + '+' + ps['name'] + '-3)/6)*' + str(prob);
                            st = '[]'+ pf['name']+' >-1 &'+pf['name']+' < 3 & ('+ ps['name']+'=1|'+ps['name']+'=2)';
                            st2 = '[]'+ pf['name']+' > 2  & ('+ ps['name']+'=3|'+ps['name']+'=4)';
                            f.write(st + ' -> ' + probstr + ':(' + ps['name'] + '\'=' + ps[
                                'name'] + '+2) +' + '1-(' + probstr + ')' + ':(' + ps['name'] + '\'=0' + ');\n');
                            f.write(st2 + ' -> ' + probstr + ':(' + ps['name'] + '\'=' + ps[
                                'name'] + '-2) +' + '1-(' + probstr + ')' + ':(' + ps['name'] + '\'=0' + ');\n');

        for n in firstsp:
            if((n not in secondsp) and (n in all_vars) and (n not in last_points.keys())):
                if(n not in groups):
                    if ( n in proteins) :
                        probstr = '((' + firstsp[n]['name'] + '+1)/6)*0.95*' + str(firstsp[n]['prob']);
                    else:
                        probstr = '((' + firstsp[n]['name'] + '+1)/6)*0.95';
                    f.write(
                        '[] ' + firstsp[n]['name'] + '=1|' + firstsp[n]['name'] + '=2 ->' + probstr + ':(' + firstsp[n][
                            'name'] + '\'=' + firstsp[n]['name'] + '+2)+' + '1-(' + probstr + ')' + ':(' + firstsp[n][
                            'name'] + '\'=0'  + ');\n');
                else:
                    cmp = groups[n]['components'];
                    fst = '';
                    sst = '';
                    sst2 = '';
                    first = True;
                    prob = 0;
                    for c in cmp :
                        if c in proteins:
                            if(first==True):
                                fst = '('+proteins[c]['name']+'=1|'+proteins[c]['name']+'=2)';
                                sst = '('+proteins[c]['name']+'\'='+proteins[c]['name']+'+2)';
                                #sst2 = '(' + proteins[c]['name'] + '\'=' + proteins[c]['name'] +')';
                                sst2 = '(' + proteins[c]['name'] + '\'=0' + ')';
                                first = False;
                            else:
                                fst += '& ('+proteins[c]['name']+'=1|'+proteins[c]['name']+'=2)';
                                sst += '& (' + proteins[c]['name']+'\'='+proteins[c]['name']+'+2)';
                                #sst2 += '& (' + proteins[c]['name'] + '\'=' + proteins[c]['name'] + ')';
                                sst2 += '& (' + proteins[c]['name'] + '\'=0'  + ')';
                            prob = prob + proteins[c]['prob'];

                    prob = prob / 6*len(cmp);
                    probstr = str(prob*0.95);

                    f.write('[] '+fst+ '->'+probstr+':'+sst+'+'+'1-('+probstr+'):'+sst2+';\n');

        first_of_actions = [];
        for k in actions.keys():
            if k in self.need_action.values():
                elm = actions[k];
                for a in elm:
                    ent1 = actions[k][a]['entry1'];
                    if(ent1 in all_vars and k in all_vars_names):
                        if(proteins.has_key(ent1)):
                            f.write('[] '+proteins[ent1]['name']+'=3|'+proteins[ent1]['name']+'=4 -> 0.9:('+k+'\'='+proteins[ent1]['name']+' )+0.1:true;\n');

                            if(ent1 not in firstsp)and(ent1 not in secondsp):
                                f.write('[]'+proteins[ent1]['name']+'=1|'+proteins[ent1]['name']+'=2 ->0.8:('+proteins[ent1]['name']+'\'='+proteins[ent1]['name']+'+2)+0.2:('+proteins[ent1]['name']+'\'=0);\n');
                        else:
                            cmps = pathway.get_compounds();
                            if(cmps.has_key(ent1)):
                                f.write('[] '+cmps[ent1]['name']+'=3|'+cmps[ent1]['name']+'=4 -> 0.9:('+k+'\'='+cmps[ent1]['name']+' )+0.1:true;\n');

        f.write('endmodule\n');
        f.close();

        return filename;


    #calculate the last points of trim points probability
    def calcute_trim_prob(self, old_last_points={}, new_last_points={}, trim_points=[], pathway=None, degs=None, all_genes=None):

        relations = pathway.get_relations();
        proteins = pathway.get_proteins();
        components = pathway.get_compounds();
        actions = pathway.get_actions();
        groups = pathway.get_groups();


        filename = self.build_prism_model(pathway, trim_points, old_last_points, degs, all_genes);

        #build properties
        for pointid in new_last_points:

            if(str(pointid) in proteins.keys()):
                vname = proteins[str(pointid)]['name'];
                pro = 'P=? [F<1000(' + vname + '= 3 )]'
                profile = 'prop.csl';
                prof = open(self.prism_path+'/'+profile, 'w');
                prof.write(pro+'\n');
                prof.close();

                res = self.run_prism_file_in_windows(filename, 'prop.csl');
                if (res != None):
                    new_last_points[pointid]['prob4'] = res;
                else:
                    new_last_points[pointid]['prob4'] = '0.0';

                pro = 'P=? [F<1000(' + vname + '= 4 )]'
                profile = 'prop.csl';
                prof = open(self.prism_path + '/' + profile, 'w');
                prof.write(pro + '\n');
                prof.close();

                res = self.run_prism_file_in_windows(filename, 'prop.csl');
                if (res != None):
                    new_last_points[pointid]['prob4'] = res;
                else:
                    new_last_points[pointid]['prob4'] = '0.0';









