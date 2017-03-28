#!/usr/bin/python

import os

with open("%s/masterBetaRunList.txt"%os.getenv("OCTET_LIST"), "wt") as out:
    with open("%s/OctetList_20112012.txt"%os.getenv("OCTET_LIST")) as fp:
	for line in fp:
	    if (line[:6]=="Octet:"):
		words=line.split()
		for idx, w in enumerate(words):
                    if w[0]=='A' or w[0]=='B':
                        out.write(words[idx+2] + ' ' + w + '\n')

    with open("%s/OctetList_20122013.txt"%os.getenv("OCTET_LIST")) as fp:
	for line in fp:
	    if (line[:6]=="Octet:"):
		words=line.split()
		for idx, w in enumerate(words):
                    if w[0]=='A' or w[0]=='B':
                        out.write(words[idx+2] + ' ' + w + '\n')
                    
