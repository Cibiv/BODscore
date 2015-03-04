#!/usr/bin/python3
import subprocess

CMD_PATH = '~/localdata/BODscore/bin/vshape-1.0.0/'
REFERENCE = '../reference/TAIR10.fas'
DATABASE = 'vcf.db'
VCF_PATH = '~/localdata/vcf/'
BAM_PATH = '~/localdata/'

base_cmd = CMD_PATH + './vshape -v -b 10 -x mitochondria -x chloroplast -l 50 -r ' + REFERENCE + ' -d ' + DATABASE 
# ' ref -q HL7.bam -s "${BASE_DIR}${List1[i]}"  -d "$DB" -t "${TAGS[i]}"  '

# ['tag': '' , 'vcf': '', 'bam': '']
samples = [
# {'tag': 'HL10_ngm_nm4' , 'vcf': 'HL10-rmdup-clipOverlap-q20-nm4-ems-annotation.vcf', 'bam': 'HL10-rmdup-clipOverlap-nm4-q20.bam'},
#{'tag': 'ABD241_bwa_nm3' , 'vcf': '20130426.B-ABD241-bwa-rmdup-clipOverlap-q20-nm3-ems-annotation.vcf', 'bam': '20130426.B-ABD241-bwa-rmdup-clipOverlap-q20-nm3.bam'},
#{'tag': 'ABD173_bwa_nm3' , 'vcf': '20130426.B-ABD173-bwa-rmdup-clipOverlap-q20-nm3-ems-annotation.vcf', 'bam': '20130426.B-ABD173-bwa-rmdup-clipOverlap-q20-nm3.bam'},
{'tag': 'ABD192_bwa_nm3' , 'vcf': '20140516.B-ABD192-bwa-rmdup-clipOverlap-q20-nm3-ems-annotation.vcf', 'bam': '20140516.B-ABD192-bwa-rmdup-clipOverlap-q20-nm3.bam'},
#{'tag': 'ABD159_bwa_nm3' , 'vcf': '20120427-ABD159-bwa-rmdup-clipOverlap-q20-nm3-ems-annotation.vcf', 'bam': '20120427-ABD159-bwa-rmdup-clipOverlap-q20-nm3.bam'},
]

def flags(d, BAM_PATH, VCF_PATH):
    return ' -q %s ' % (BAM_PATH + d['bam']) + ' -s %s ' % (VCF_PATH + d['vcf']) + ' -t %s ' % d['tag']

for ss in samples:
    cmd = base_cmd + flags(ss, BAM_PATH, VCF_PATH)
    print(cmd)
    p = subprocess.Popen(cmd, shell=True)
    p.wait()
    print("-------------------------------------------------->")
