__author__ = 'Mansoori'

from os import listdir
from os.path import isfile, join


result_path = './data/'
temp_random_path = './data/temp'
random_path ='./data/null_dist'
result_file = 'result'
random_file ='random'



class sorting():

    def sort(self,out_file):
        results = {}
        result_raw = [line.rstrip('\n') for line in open(join(result_path, result_file), 'r')]
        randomfiles = [f for f in listdir(random_path) if isfile(join(temp_random_path, f))]
        rlen = len(result_raw);
        cnt = 0;
        main_probs = {};
        while(cnt<rlen):
            splitr = result_raw[cnt].split('\t');
            main_probs[splitr[0]] = (float(splitr[2]));
            cnt += 1

        for fr in randomfiles:
            random_raw = [line.rstrip('\n') for line in open(join(temp_random_path, f), 'r')]
            rlen = len(random_raw);
            splitr = random_raw[0].split('\t');
            title = splitr[0] + '\t' + splitr[1];
            idp = splitr[0]
            cnt = 0;
            cntp = 0;
            while(cnt < rlen):
                splitr = random_raw[cnt].split('\t');
                prob = float(splitr[0])
                if(prob > main_probs[idp]):
                    cntp+=1;
                cnt+=1;
            pvalue = float(cntp) / rlen;
            if(results.has_key(pvalue)==False):
                results[pvalue] = []
            results[pvalue].append({'title':title, 'score':main_probs[idp], 'pvalue':float(cntp)/rlen});

        #sort based on the pvalue and score
        pkeys = results.keys();
        for i in range(0, len(pkeys)):
            for j in range(i, len(pkeys)):
                if(pkeys[i]>pkeys[j]):
                    temp = pkeys[i];
                    pkeys[i] = pkeys[j];
                    pkeys[j] = temp;

        for pk in pkeys:
            pkpvalues = results[pk];
            for i in range(0, len(pkpvalues)):
                for j in range(i, len(pkpvalues)):
                    if(pkpvalues[i]['score'] < pkpvalues[j]['score']):
                        temp = pkpvalues[i];
                        pkpvalues[i] = pkpvalues[j];
                        pkpvalues[j] = temp;

        fout = open(join(result_path, out_file), 'w')
        for pk in pkeys :
            pkpvalues = results[pk];
            for pkp in pkpvalues:
                fout.write(pkp['title']+'\t'+str(pkp['score'])+'\t'+str(pkp['pvalue'])+'\n');


    def sort_stored_random(self, out_file):
        results = {}
        result_raw = [line.rstrip('\n') for line in open(join(result_path, result_file), 'r')]
        randomfiles = [f for f in listdir(random_path) if isfile(join(random_path, f))]
        rlen = len(result_raw);
        cnt = 0;
        main_probs = {};
        while (cnt < rlen):
            splitr = result_raw[cnt].split('\t');
            main_probs[splitr[0]] = (float(splitr[2]));
            cnt += 1

        for fr in randomfiles:
            random_raw = [line.rstrip('\n') for line in open(join(random_path, fr), 'r')]
            rlen = len(random_raw);
            splitr = random_raw[0].split('\t');
            title = splitr[0] + '\t' + splitr[1];
            idp = splitr[0]
            cnt = 0;
            cntp = 0;
            while (cnt < rlen):
                splitr = random_raw[cnt].split('\t');
                prob = float(splitr[2])
                mprob = 0;
                if (idp in main_probs):
                    mprob = main_probs[idp];
                if (prob >= mprob):
                    cntp += 1;
                cnt += 1;
            pvalue = float(cntp) / rlen;
            if (results.has_key(pvalue) == False):
                results[pvalue] = []
            results[pvalue].append({'title': title, 'score': mprob, 'pvalue': pvalue});

        # sort based on the pvalue and score
        pkeys = results.keys();
        for i in range(0, len(pkeys)):
            for j in range(i, len(pkeys)):
                if (pkeys[i] > pkeys[j]):
                    temp = pkeys[i];
                    pkeys[i] = pkeys[j];
                    pkeys[j] = temp;

        for pk in pkeys:
            pkpvalues = results[pk];
            for i in range(0, len(pkpvalues)):
                for j in range(i, len(pkpvalues)):
                    if (pkpvalues[i]['score'] < pkpvalues[j]['score']):
                        temp = pkpvalues[i];
                        pkpvalues[i] = pkpvalues[j];
                        pkpvalues[j] = temp;

        fout = open(join(result_path, out_file), 'w')
        for pk in pkeys:
            pkpvalues = results[pk];
            for pkp in pkpvalues:
                fout.write(pkp['title'] + '\t' + str(pkp['score']) + '\t' + str(pkp['pvalue']) + '\n');

