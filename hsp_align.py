from Bio.Blast import NCBIXML
from Bio import Seq, SeqRecord, SeqIO
import csv, sys, os

def hsp_align(blast_xml, outfile, hub_strain_fasta):
    result_handle = open(blast_xml)
    blast_records = NCBIXML.parse(result_handle)

    e_val_threshold = 0.00001

    with open(outfile, 'w') as outfile:
        fieldnames = ['hsps', 'query_id', 'subj_accession', 'subj_name', 'query_len', 'pident', 'e-val', 'subj_start', 'subj_end', 'subj_seq']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                alignment_dict = dict()
                hsps_meeting_threshold = []
                hsp_ct = 0
                hsp_coords = []
                for hsp in alignment.hsps:
                    if hsp.expect < e_val_threshold:
                        hsp_ct += 1
                        hsps_meeting_threshold.append(hsp)
                        hsp_coords.append(hsp.sbjct_start)
                        hsp_coords.append(hsp.sbjct_end)
                if hsp_ct > 0:
                    alignment_dict['subj_seq'] = [hsp.sbjct for hsp in hsps_meeting_threshold][0]
                    alignment_dict['hsps'] = hsp_ct
                    alignment_dict['query_id'] = blast_record.query
                    alignment_dict['subj_accession'] = alignment.hit_id.split('|')[1]
                    alignment_dict['subj_name'] = alignment.hit_def.split(' ')[0][0] + '_' + alignment.hit_def.split(' ')[1]
                    alignment_dict['query_len'] = blast_record.query_length
                    alignment_dict['e-val'] = [hsp.expect for hsp in hsps_meeting_threshold][0]
                    alignment_dict['subj_start'] = min(hsp_coords)
                    alignment_dict['subj_end'] = max(hsp_coords)
                    for good_hsp in hsps_meeting_threshold:
                        mismatch = good_hsp.align_length - good_hsp.identities - good_hsp.gaps
                        pident = (good_hsp.align_length - mismatch - good_hsp.gaps)/good_hsp.align_length
                    alignment_dict['pident'] = pident * 100
                if hsp_ct > 1:
                    alignment_dict['e-val'] = [hsp.expect for hsp in hsps_meeting_threshold]
                    alignment_dict['pident'] = None
                    with open(hub_strain_fasta, 'r') as fasta_in:
                        record_dict = SeqIO.to_dict(SeqIO.parse(fasta_in, 'fasta'))
                        subj_seq = record_dict[alignment_dict['subj_accession']].seq[(alignment_dict['subj_start']-1):alignment_dict['subj_end']]
                        alignment_dict['subj_seq'] = subj_seq
                if hsp_ct > 0:
                    writer.writerow(alignment_dict)
                alignment_dict.clear()

if __name__ == '__main__':
    if len(sys.argv) == 4:
         hsp_align(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
         print("Usage: hsp_align.py blast_in.xml alignment_out.csv hub_strains.fa")
         sys.exit(0)
