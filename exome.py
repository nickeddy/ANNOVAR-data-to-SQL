#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, csv, re, subprocess
from optparse import OptionParser

version = "%prog v.1"
description = "This script creates a database of exomic DNA data via ANNOVAR files."
rs_regex = "(?<=Name\=)rs\d+"
rs_regex_obj = re.compile(rs_regex)

parser = OptionParser(version=version, description=description)
parser.add_option("--path", dest="path", help="The path to the folder containing all of the ANNOVARin files. Needs trailing slash.", metavar="</path/to/files/>")
parser.add_option("--name", dest="name", help="The name preceding the .ANNOVARin in the filenames e.g. raw_variants", metavar="<common_name>")
parser.add_option("--output", dest="output", help="The output file for SQL query", metavar="</path/to/output>")
(options, args) = parser.parse_args()

mandatory_options = ['path', 'name', 'output']
for m in mandatory_options:
    if not options.__dict__[m]:
        parser.print_help()
        exit(-1)

if not os.path.exists(options.__dict__['path']):
    print "\n\033[1;31m%s does not exist!\033[1;m\n" % (options.__dict__['path'])
    parser.print_help()
    exit(-1)

if not os.path.exists(options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.variant_function'):
    print "\n\033[1;31m%s does not exist!\033[1;m\n" % (options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.variant_function')
    parser.print_help()
    exit(-1)

variant_function = open(options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.variant_function')
variant_reader = csv.reader(variant_function, delimiter='    ')

items = dict()
for row in variant_reader:
    items[row[12]] = [row[2], row[3], row[4], row[9], row[8], row[10], row[11], row[7], row[0], row[0], row[5], row[6]]

if not os.path.exists(options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.hg19_snp135'):
    print "\n\033[1;31m%s does not exist!\033[1;m\n" % (options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.hg19_snp135')
    parser.print_help()
    exit(-1)

# http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=rs76987352
hg19_snp135 = open(options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.hg19_snp135')
hg19_reader = csv.reader(hg19_snp135, delimiter='    ')

for row in hg19_reader:
    rs_match = rs_regex_obj.findall(row[1])
    items[row[12]].append(rs_match[0])

for item in items:
    try:
        test = items[item][12]
    except:
        items[item].append("N/A")
if not os.path.exists(options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.hg19_ALL.sites.2010_11_dropped'):
    print "\n\033[1;31m%s does not exist!\033[1;m\n" % (options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.hg19_ALL.sites.2010_11_dropped')
    parser.print_help()
    exit(-1)

_2010_11_dropped = open(options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.hg19_ALL.sites.2010_11_dropped')
_2010_11_reader = csv.reader(_2010_11_dropped, delimiter='    ')

for row in _2010_11_reader:
    items[row[12]].append(row[1])

for item in items:
    try:
        test = items[item][13]
    except:
        items[item].append("N/A")

if not os.path.exists(options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.hg19_ALL.sites.2011_05_dropped'):
    print "\n\033[1;31m%s does not exist!\033[1;m\n" % (options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.hg19_ALL.sites.2011_05_dropped')
    parser.print_help()
    exit(-1)

_2011_05_dropped = open(options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.hg19_ALL.sites.2011_05_dropped')
_2011_05_reader = csv.reader(_2011_05_dropped, delimiter='    ')

for row in _2011_05_reader:
    items[row[12]].append(row[1])

for item in items:
    try:
        test = items[item][14]
    except:
        items[item].append("N/A")

if not os.path.exists(options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.exonic_variant_function'):
    print "\n\033[1;31m%s does not exist!\033[1;m\n" % (options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.exonic_variant_function')
    parser.print_help()
    exit(-1)

evariant = open(options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.exonic_variant_function')
evariant_reader = csv.reader(evariant, delimiter='    ')

for row in evariant_reader:
    items[row[13]].append(row[2])
    items[row[13]].append(row[1])

for item in items:
    try:
        test = items[item][15]
    except:
        items[item].append("N/A")
for item in items:
    try:
        test = items[item][16]
    except:
        items[item].append("N/A")

hg19_ljb_pp2 = open(options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.hg19_ljb_pp2_dropped')
hg19_ljb_pp2_reader = csv.reader(hg19_ljb_pp2, delimiter='    ')

for row in hg19_ljb_pp2_reader:
    items[row[12]].append(row[1])

for item in items:
    try:
        test = items[item][17]
    except:
        items[item].append("N/A")

hg19_avsift = open(options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.hg19_avsift_dropped')
hg19_avsift_reader = csv.reader(hg19_avsift, delimiter='    ')

for row in hg19_avsift_reader:
    items[row[12]].append(row[1])

for item in items:
    try:
        test = items[item][18]
    except:
        items[item].append("N/A")

hg19_ljb_mt = open(options.__dict__['path'] + options.__dict__['name'] + '.ANNOVARin.hg19_ljb_mt_dropped')
hg19_ljb_mt_reader = csv.reader(hg19_ljb_mt, delimiter='    ')

for row in hg19_ljb_mt_reader:
    items[row[12]].append(row[1])

for item in items:
    try:
        test = items[item][19]
    except:
        items[item].append("N/A")


# the following code is probably pretty terrible.
sql = "CREATE TABLE " + options.__dict__['name'] + " (id int NOT NULL AUTO_INCREMENT, chr varchar(6), start int, end int, coverage smallint, var_conf_score decimal(7,2), map_quality decimal(4,2), quality_by_depth decimal(4,2), state varchar(3), gene_type varchar(20), exon_intron varchar(20), ref_base varchar(1), alt_base varchar(1), rs_num varchar(15), minor_allele_freq decimal(10,10), minor_allele_freq_source decimal(10,10), protein_change blob, name_of_change varchar(23), polyphen_score decimal(10,10), sift_score decimal(10,10), mutation_taster decimal(10,10), PRIMARY KEY(id));\n"
for item in items:
    _1 = "N/A"
    _2 = "N/A"
    _3 = "N/A"
    _4 = "N/A"
    _5 = "N/A"
    _6 = "N/A"
    _7 = "N/A"
    _8 = "N/A"
    _9 = "N/A"
    _10 = "N/A"
    _11 = "N/A"
    _12 = "N/A"
    _13 = "N/A"
    _14 = "N/A"
    _15 = "N/A"
    _16 = "N/A"
    _17 = "N/A"
    _18 = "N/A"
    _19 = "N/A"
    _20 = "N/A"
    try:
        _1 = items[item][0]
    except:
        pass
    try:
        _2 = items[item][1]
    except:
        pass
    try:
        _3 = items[item][2]
    except:
        pass
    try:
        _4 = items[item][3]
    except:
        pass
    try:
        _5 = float(items[item][4])
    except:
        pass
    try:
        _6 = float(items[item][5])
    except:
        pass
    try:
        _7 = float(items[item][6])
    except:
        pass
    try:
        _8 = items[item][7]
    except:
        pass
    try:
        _9 = items[item][8]
    except:
        pass
    try:
        _10 = items[item][9]
    except:
        pass
    try:
        _11 = items[item][10]
    except:
        pass
    try:
        _12 = items[item][11]
    except:
        pass
    try:
        _13 = items[item][12]
    except:
        pass
    try:
        _14 = float(items[item][13])
    except:
        pass
    try:
        _15 = float(items[item][14])
    except:
        pass
    try:
        _16 = items[item][15]
    except:
        pass
    try:
        if re.match("\d", items[item][16]):
            _17 = "N/A"
        else:
            _17 = items[item][16]
    except:
        pass
    try:
        _18 = float(items[item][17])
    except:
        pass
    try:
        _19 = float(items[item][18])
    except:
        pass
    try:
        _20 = float(items[item][19])
    except:
        pass
    sql += "INSERT INTO %s (chr, start, end, coverage, var_conf_score, map_quality, quality_by_depth, state, gene_type, exon_intron, ref_base, alt_base, rs_num, minor_allele_freq, minor_allele_freq_source, protein_change, name_of_change, polyphen_score, sift_score, mutation_taster) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s');\n" % (options.__dict__['name'], _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20)

#	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----
#		0		1		2		3		4		5		6		7		8		9		10		11		12		13		14		15		16		17		18		19
#	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----
#		chr		start	end		cov		v_c_s	map_q	q_b_d	state	gene	exon	ref_b	alt_b	rs_num	maf		mafs	p_c		n_o_c	pp_s	s_s		m_t
#	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----	----
with open(options.__dict__['output'], 'w') as output_file:
    output_file.write(sql)
