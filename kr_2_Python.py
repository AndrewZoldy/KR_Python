from Bio import SeqIO
import re
from Bio.Blast import NCBIWWW
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
import argparse
import os
import time
import urllib
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='auto blast name giver')
    parser.add_argument("-i", "--input_file", type=str, help='path to input file', default=os.curdir+'/pacbio-only-smrt.fasta')
    parser.add_argument('-fr', '--fresh_start', action='store_true', help='remove')

    args = parser.parse_args()
    inp = args.input_file
    fresh = args.fresh_start
    spec_finder = re.compile(r'''\|\s(PREDICTED:\s)*(?P<name>[\w]*\s[\w]*)[\s\n]''')
    new_list_seqs = []
    record_list = []
    count_list = []
    real_list_bro = []
    if fresh == True:
        os.remove('session.txt')
    else:
        pass
    k = 0
    with open(inp, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            record_list.append(record)
        try:
            with open('session.txt', 'r') as f:
                k = json.load(f)
            print('Search will start with the last work progress')
        except json.decoder.JSONDecodeError:
            print('Start writing the dump.file')
        except FileNotFoundError:
            print('Initiation of the dump.file')
        while k <len(record_list):
            try:
                record = record_list[k]
                print('Searching in NCBI blast...')
                res = NCBIWWW.qblast('blastn', 'nt', record.seq,hitlist_size=5,alignments=5,descriptions=5,
                                     format_type='XML')
                print('QBLASTed successfully')
                par = NCBIXML.parse(res)
                # print('NAME___',record.id)
                for i in par:
                    for al in i.alignments:
                        hsp_string=''
                        exp_list=[]
                        for hsp in al.hsps:
                            exp_list.append(hsp.expect)
                        exp=min(exp_list)
                        founded = spec_finder.finditer(al.title)
                        #print(al.title)
                        for j in founded:
                            record_features = [exp,j.group('name')]
                            new_list_seqs.append(record_features)
                new_list_seqs=sorted(new_list_seqs)
                temp_list=[]
                temp_list.append(new_list_seqs[0])
                for i in range(len(new_list_seqs)-1):
                    if (new_list_seqs[i+1][0]-new_list_seqs[i][0])>0.000001:
                        break
                    else:
                        temp_list.append(new_list_seqs[i+1])
                for t in temp_list:
                    count_list.append(t[1])
                names_set=set(count_list)
                names_set=list(names_set)

                for i in range(len(names_set)):
                    x = count_list.count(names_set[i])
                    names_set[i] = [x,names_set[i]]
                name = max(names_set)[1]
                print('Species of this contig is:',name)
                named_record = SeqRecord(record.seq, name, '', '')
                real_list_bro.append(named_record)
                new_list_seqs = []
                count_list = []
                k+=1
                with open("session.txt", "w") as f:
                    json.dump((k), f)

                time.sleep(5)

                SeqIO.write(real_list_bro, inp.replace('.fasta', '') + '_temp_out.fasta', 'fasta')
            except urllib.error.URLError:
                print('connection error, wait please')
                time.sleep(5)
                continue
            except ValueError:
                print('connection error, wait please')
                time.sleep(5)
                continue
    names=[]
    with open(inp.replace('.fasta', '') + '_temp_out.fasta','r') as handle:
        data=SeqIO.parse(handle, 'fasta')
        for record in data:
            names.append(record.id)
        names=set(names)
        names=list(names)
        payload=[]
        for name in names:
            for record in data:
                if record.id == name:
                    payload.append(record)
        SeqIO.write(payload, inp.replace('.fasta', '') + '_out.fasta')
    os.remove(inp.replace('.fasta', '') + '_temp_out.fasta')
    os.remove('session.txt')
