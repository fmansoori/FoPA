__author__ = 'Mansoori'

import random





class prism_trim():
    rgen = random.Random();

    prism_path = '..\data\prism output\p1var';

    m = 7


    first_of_actions = []
    seconds = {};
    firsts = {};

    first_of_pathways = []
    all_vars = [];

    not_analyzed_last_points = {}
    analyzed_last_points = {}
    analyzed_vars = []

    def __init__(self):
        self.first_of_actions = []
        self.seconds = {};
        self.firsts = {};

        self.first_of_pathways = []
        self.all_vars = [];

        self.not_analyzed_last_points = {}
        self.analyzed_last_points = {}
        self.analyzed_vars = []
        self.m = 7;


    def find_cycle(self, root, traverse=[]):

        traverse.append(root);
        root_relations = [];
        if(self.firsts.has_key(root)):
            root_relations = self.firsts[root];
        for rr in root_relations:
            sr = root_relations[rr]['second'];
            if(sr in traverse):
                return True;
            elif (sr in self.seconds):
                if(self.find_cycle(sr) == True):
                    return True;

        return False;

    def find_path(self, root, fvars=[], seconds={}, traverse=[]):

        skeys = seconds.keys();
        if (root in skeys):
            root_relations = seconds[root];
            for fri in root_relations:
                fr = root_relations[fri]['first'];
                if (fr in skeys):
                    if(fr not in traverse):
                        traverse.append(root);
                        self.find_path(fr, fvars, seconds, traverse);
                else:
                    if(fr not in fvars):
                        fvars.append(fr);
        else:
            if(root not in fvars):
                fvars.append(root);


    def find_firsts(self, pfirsts):
        max = 0;
        #find the node with maximum input
        for s in self.seconds:
            slen = len(self.seconds[s])
            if(max < slen):
                mid = s;
                max = slen;

        if(len(self.seconds)!= 0):
            self.find_path(mid, pfirsts, self.seconds);
            if (max+2) > 10:
                self.m = max+2;
            else:
                self.m = 6;

    def set_all_variables(self, vars=[]):

        self.all_vars = vars;

    def set_firsts_seconds(self, firsts={}, seconds={}):

        self.firsts = firsts.copy();
        self.seconds = seconds.copy();

        self.cal_first_of_pathway();


    def cal_first_of_pathway(self):
        self.first_of_pathways = [];
        for f in self.firsts.keys():
            if f not in self.seconds.keys():
                self.first_of_pathways.append(f);

    #remove relations of begins
    def remove_rel(self, begin):

        if(self.seconds.has_key(begin)):
            ss = self.seconds[begin];
            for s in ss:
                fi = ss[s]['first'];

                if(fi in self.firsts):
                    ff = self.firsts[fi];
                    ks = ff.keys();
                    for cnt in ks:
                        if(ff[cnt]['second']==begin):
                            self.firsts[fi].pop(cnt);

                    if(self.firsts[fi]) == {} :
                        self.firsts.pop(fi);

            self.seconds.pop(begin);

    #check to see if the second of the begin points are in the trim points if not remove the begin from trim points
    def trim_investigate(self, begin, trim_points=[]):

        ff = self.firsts[begin];
        ks = ff.keys();
        for cnt in ks:
            if (ff[cnt]['second'] in trim_points):
                return True;

        return False;

    #def cal_trim(self, vars=[], rrelations=[], seconds=[], trim_points=[], last_points={} ):
    def trim(self, old_last_points={}, begins = {}, trim_points=[], new_last_points={} ):

        for pidx in old_last_points:
            self.not_analyzed_last_points[pidx]={'prob0:':old_last_points[pidx]['prob0'], 'prob1': old_last_points[pidx]['prob1'],\
                'prob2':old_last_points[pidx]['prob2'], 'prob3':old_last_points[pidx]['prob3'], 'prob4':old_last_points[pidx]['prob4']};

        ppoints = []

        #potansel first points
        pfirsts1 = [];
        self.find_firsts(pfirsts1);

        for a in pfirsts1:
            if (a not in self.all_vars) or (a in self.analyzed_vars):
                pfirsts1.pop(a);

        pfirsts2 = []
        self.cal_first_of_pathway();
        for a in self.first_of_pathways:
            if (a in self.all_vars) and (a not in self.analyzed_vars) and (a not in pfirsts1):
                pfirsts2.append(a);

        if (len(pfirsts1) > 0 or len(pfirsts2) > 0) :
            m = self.m;
            ti = False;
            while(len(trim_points)<=1 or ti==False):
                if(len(pfirsts1) != 0):
                    #select one element of pfirsts randomly
                    #rf = self.rgen.randint(0, len(pfirsts) - 1);
                    b = pfirsts1.pop();

                elif(len(pfirsts2) != 0):
                    b = pfirsts2.pop();
                else:
                    break;

                new_last_points.clear();
                stack = [];
                self.cal_trim(b, m, trim_points, ppoints, new_last_points, stack);

                ti = self.trim_investigate(b, trim_points);
                if(ti == False):
                    m = self.m;
                    while(len(ppoints)!=0):
                        ppoints.pop();
                    while(len(trim_points)!=0):
                        trim_points.pop();

                if(len(trim_points)==1)and(len(pfirsts1)==0and(len(pfirsts2)==0)):
                    trim_points.pop();

                print '.',

            if b in self.not_analyzed_last_points:
                begins[b] = self.not_analyzed_last_points[b];
            else:
                begins[b] = {'prob0':'0.0', 'prob1':'0.0', 'prob2':'0.0', 'prob3':'0.0', 'prob4':'0.0'}

            for t in trim_points:
                self.remove_rel(t);
                #if(t not in new_last_points.keys()):
                #    self.analyzed_vars.append(t);
                if (t in self.not_analyzed_last_points.keys()) and (t not in new_last_points.keys()):
                    self.analyzed_last_points[t] = self.not_analyzed_last_points[t];
                    self.not_analyzed_last_points.pop(t);
                if t in self.first_of_pathways:
                    if t in self.not_analyzed_last_points:
                        begins[t]= self.not_analyzed_last_points[t];
                    else:
                        begins[t] = {'prob0':'0.0', 'prob1':'0.0', 'prob2':'0.0', 'prob3':'0.0', 'prob4':'0.0'};
                if(t in self.seconds ) and (t not in self.analyzed_last_points):
                    sl = self.seconds[t];
                    for elm in sl:
                        s = sl[elm]['first'];
                        if (s not in trim_points) and (s in self.analyzed_last_points):
                            begins[s] = self.analyzed_last_points[s];
                            trim_points.append(s);



    def cal_trim(self, begin, d, trim_points=[], ppoints=[], last_points={}, stack=[] ):

        d = min(d, self.m-len(trim_points));
        if (d <= 0) :
            return False;

        sb = {}
        if (begin in self.seconds.keys()):
            sb = self.seconds[begin];

        can = True;

    #    td = d-1;
        ppoints.append(begin);
        stack.append(begin);
        cnt = 0;
        while(can != False) and (cnt < len(sb) and d > 0):
            ssb = sb[cnt];
            if (ssb['first'] in self.all_vars) and (ssb['first'] not in self.analyzed_vars):
                if (ssb['first'] not in trim_points):
                #if(ssb['first'] not in trim_points) and ssb['first'] not in ppoints:
                    #check for cycle
                    if(ssb['first'] not in stack) :
                        if(self.cal_trim(ssb['first'], d, trim_points, ppoints, last_points) == False):
                            can = False;
                        else :
                            d -= 1;
                            d = min(d, self.m-len(trim_points));

            print '.',
            cnt+=1;

        if(can == True and cnt == len(sb)):
            trim_points.append(begin);
            if(begin in ppoints):
                ppoints.remove(begin);
            d-=1;


        stack.remove(begin);
        d = min(d, self.m-len(trim_points));

        had = False;
        if(len(ppoints)==0):
            fa = {}
            if(d > 0 ):
                if(begin in self.firsts):
                    fa = self.firsts[begin];

                for cnt in fa.keys():
                    if(d> 0):
                        ffa = fa[cnt]
                        if ((ffa['second'] in self.all_vars)) and (ffa['second'] not in trim_points) and (ffa['second']not in self.analyzed_vars):
                            if ffa['second'] not in ppoints :
                                if self.cal_trim(ffa['second'], d, trim_points, ppoints, last_points) == True:
                                    had = True;
                                    d-=1;
                                    d = min(d, self.m-len(trim_points));


        if(had == False) and (begin in trim_points):
            last_points[begin]={'prob0':'0.0','prob1':'0.0','prob2':'0.0','prob3':'0.0','prob4':'0.0'};

        if (begin in trim_points) and (begin in self.seconds):
            sa = self.seconds[begin];
            for idx in sa:
                sap = sa[idx];
                if(sap['first'] in last_points.keys()):
                    last_points.pop(sap['first']);

        if(begin in trim_points):
            return True;
        else:
            return False;
















