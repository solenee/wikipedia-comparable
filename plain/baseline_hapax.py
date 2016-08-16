#coding:utf-8

#-------------------------------
# Parse ts.xml file to save the
# test set into a CVS file
# separatation : " : "
#-------------------------------

import numpy as np
from numpy.matlib import zeros
import json
import re, sys, os
import xml.etree.ElementTree as ET
from codecs import open
import time
from collections import defaultdict, Counter

HAPAX = {} #inverted index
FR="src"
EN="trgt"
IDs = { FR : {}, EN : {} } #annuary of src, tgrt files
currentID = { FR : -1 , EN : -1 }
INDEX = {}
INDEX_EN = defaultdict(list)

S_TOP = "top"
S_PER = "percentage"
SELECTION = S_PER
TOP = 1000
MAX_TOP = 1000
PER = 0.9

RES_FILE = "baseline.res90p"


STATE = None #matching matrix 
nbIter = 0
eps = 0.0000001 #np.finfo(float).eps #precision
alpha = 1 #0000

def loadJson(filename) : 
  with open(filename, "r", encoding="utf-8") as f:
    text = f.read()
    jsonDesc = json.loads(text)
  return jsonDesc

def savestate() :
  with open('bckp/invertedIndex.txt', 'w') as outfile:
      json.dump(HAPAX, outfile)
      print "Bckp inverted index state in bckp/invertedIndex.txt"
  with open('bckp/annuary.txt', 'w') as outfile:
    json.dump(IDs, outfile)
    print "Bckp IDs saved bckp/annuary.txt"
      
def readFiles(directory, language) :
  """" language : src or trgt. Index the files contained in directory """
  start_time = time.time()
  file_paths = os.listdir(directory)
  if language == EN :
    fr_hapaxes = set([item for sublist in loadJson('indexFR.txt').values() for item in sublist]) #reduce(lambda a, b : set(a) | set(b), loadJson('indexFR.txt').values(), [])
    print len(fr_hapaxes)
  for file_path in file_paths:
    #Give an ID
    currentID[language] = currentID[language] + 1
    IDs[language][currentID[language]] = file_path
    fileID = currentID[language]
    #Index file
    if language == FR : INDEX[file_path] = index(directory+"/"+file_path, language)
    elif language == EN : 
      hapax_list = index(directory+"/"+file_path, language)
      #Update inverted index
      for h in hapax_list :
        if h in fr_hapaxes : INDEX_EN[h].append(fileID)
##        if h not in HAPAX.keys() :
##          HAPAX[h] = {}
##          HAPAX[h][FR] = []
##          HAPAX[h][EN] = []
##        HAPAX[h][language].append(file_path)
##      elapsed_time = time.time() - start_time
##      if (elapsed_time > 3600 ) :
##        print "Saving state"
##        #savestate()
##        print "Current state saved"
##        start_time = time.time()
  elapsed_time = time.time() - start_time
  print str(elapsed_time)

def index(file_path, language) :
  dico = defaultdict(int)
  with open(file_path, "r", encoding="utf-8") as f: #codecs.open(file, "r", "utf-8") as f:
    tmp = re.split('(\s+)', f.read())
    for w in tmp :
      token = w.strip()
      if len(token) >= 4 : dico[token] += 1
  hapax = [w for w in dico.keys() if dico[w] == 1]
  #print hapax
  return hapax

def lookupEN(hapax_list) :
  """ Find all IDS of EN files which contain at least one item of the given list in their index """
  result = []
  for h in hapax_list :
    result.extend(INDEX_EN.get(h, []))
  scores = Counter(result)
  #print scores
  return scores


def keepCommonHapaxesOnly() :
  for hpx in HAPAX.keys() :
    if len(HAPAX[hpx][FR]) == 0 or len(HAPAX[hpx][EN]) == 0 : del HAPAX[hpx]

  
if __name__ == "__main__":
  shouldIndex = False
  if shouldIndex :
##    print "Indexing FR ..."
##    readFiles("Donnees/FR/fr", FR)
##    with open('invertedIndexFr.txt', 'w') as outfile:
##      json.dump(HAPAX, outfile)
##      print "Current state saved in invertedIndexFr.txt"
##    with open('indexFR.txt', 'w') as outfile:
##      json.dump(INDEX, outfile)
##      print "Current state saved in indexFR.txt"
##    with open('annuaryFR.txt', 'w') as outfile:
##      json.dump(IDs[FR], outfile)
##      print "IDs saved annuaryFR.txt"
##    with open('indexFR.csv', 'w', encoding="utf-8") as csvfile:
##      for idSrc in IDs[FR].keys() :
##        csvfile.write(str(idSrc)+"\t"+" ; ".join(INDEX[IDs[FR][idSrc]]))
##      print "Index FR saved in indexFR.csv"
  
    print "Indexing EN ..."
    readFiles("Donnees/EN/en", EN)
    with open('invertedIndexEn.txt', 'w', encoding="utf-8") as outfile:
      json.dump(HAPAX, outfile)
      print "Current state saved in invertedIndexEn.txt"
    with open('indexEN.txt', 'w', encoding="utf-8") as outfile:
      json.dump(INDEX_EN, outfile)
      print "Current state saved in indexEN.txt"
    with open('annuaryEN.txt', 'w', encoding="utf-8") as outfile:
      json.dump(IDs[EN], outfile)
      print "IDs saved annuaryEN.txt"
  else :
    start_time = time.time()
##    print "Loading IDs (annuary.txt) ..."
##    IDs[FR] = loadJson('annuaryFR.txt')
##    print "Nb files in FR : "+ str(len(IDs[FR].keys()))    
##    print "Loading inverted index FR(indexFR.txt) ..."
##    INDEX = loadJson('indexFR.txt')
##    print len(INDEX.keys())
##    print "\n"
##    with open('indexFR_names.csv', 'w', encoding="utf-8") as csvfile:
##      for idSrc in IDs[FR].keys() :
##        csvfile.write(IDs[FR][idSrc]+"\t"+" ; ".join(INDEX[IDs[FR][idSrc]])+"\n")
##    print "Index FR saved in indexFR_names.csv"
    
    print "Loading inverted index EN(indexEN.txt) ..."
    INDEX_EN = loadJson('indexEN.txt')
    print len(INDEX_EN.keys())
    print "\n"
    IDs[EN] = loadJson('annuaryEN.txt')
    print "Nb files in EN : "+ str(len(IDs[EN].keys()))
    elapsed_time = time.time() - start_time
    print str(elapsed_time)+"\n"

  #initialize matching matrix
  start_time = time.time()
  nbIter = 0
  nr = 114789 #len(loadJson('indexFR.txt').keys())
  nc = len(IDs[EN].keys())
  #STATE = np.zeros((nr, nc))
  thres = float(1) / (nc)
  print "threshold = "+str(thres)
  #PROB = {"0" : np.zeros((nr, nc)), "1" : np.zeros((nr, nc)) }
  #print STATE
  #print PROB
  #maxlen = 0
  trgtFiles = set()
  print "Matrix size : "+str(np.size(STATE) )
  with open(RES_FILE, 'w', encoding="utf-8") as resFile:
    with open('indexFR_names.csv', 'r', encoding="utf-8") as csvfile:
      for lineFR in csvfile :
        tmp = lineFR.split("\t")
        idSrc = tmp[0]
        hapax_list = tmp[1].split(" ; ")
  #for idSrc in IDs[FR].keys() :
    #hapax_list = index(directory+"/"+IDs[FR][idSrc], FR)
        candidates = lookupEN(hapax_list)
        #maxlen = max(maxlen, len(candidates))
        resFile.write(idSrc+"\t")#resFile.write(IDs[FR][idSrc]+"\t")
        scores = []
        ###max([targetFreqCounts[x] for x in bilingualDico[word2] if targetFreqCounts[x] > 0])
        if SELECTION == S_TOP :
          for idTrgt, count in candidates.most_common(TOP) :
            #if candidates[idTrgt] > 1 :
            scores.append(IDs[EN][str(idTrgt)]+":"+str(count))
        elif SELECTION == S_PER :
          maxx = candidates.most_common(1)
          if len(maxx) == 0 : print "no candidates found for "+str(idSrc)
          else :
            idMax, maxval = maxx[0]
            thres = maxval*PER
            #print thres
            for idTrgt, count in candidates.most_common(MAX_TOP) :
              if count >= thres :
                scores.append(IDs[EN][str(idTrgt)]+":"+str(count))
              else : break
        else :
          for idTrgt, count in candidates.most_common(MAX_TOP) :
            scores.append(IDs[EN][str(idTrgt)]+":"+str(count))
        resFile.write("\t".join(scores)+"\n")
        #print scores
  #print maxlen
        #print score
        #STATE[idSrc, idTrgt] = float(score)
        
  filename = "state_"+str(nbIter)
  print filename
  elapsed_time = time.time() - start_time
  print str(elapsed_time)+"\n"
  #print STATE
  #np.save(filename, STATE)


def next() :
  if True :
    print "Computing haxapes in common..."
    print "Total numbers of hapaxes : "+str(len(HAPAX.keys()))
    with open('RawInvertedIndex.txt', 'w') as outfile:
      json.dump(HAPAX, outfile)
      print "All happaxes saved in RawInvertedIndex.txt"
    keepCommonHapaxesOnly()
    print "Total numbers of hapaxes in common : "+str(len(HAPAX.keys()))
    #save inverted index
    with open('invertedIndex.txt', 'w') as outfile:
      json.dump(HAPAX, outfile)
      print "Final inverted index state in invertedIndex.txt"
    with open('annuary.txt', 'w') as outfile:
      json.dump(IDs, outfile)
      print "All IDs saved annuary.txt"
  else :
    start_time = time.time()
    print "Loading inverted index (invertedIndex.txt) ..."
    HAPAX = loadJson('invertedIndex.txt')
    print len(HAPAX.keys())
    print "\n"
    print "Loading IDs (annuary.txt) ..."
    IDs = loadJson('annuary.txt')
    print len(IDs.keys())
    elapsed_time = time.time() - start_time
    print str(elapsed_time)+"\n"

  #initialize matching matrix
  nbIter = 0
  nr = len(IDs[FR].keys())
  nc = len(IDs[EN].keys())
  STATE = np.zeros((nr, nc))
  thres = float(1) / (nc)
  print "threshold = "+str(thres)
  PROB = {"0" : np.zeros((nr, nc)), "1" : np.zeros((nr, nc)) }
  #print STATE
  #print PROB
  print "Matrix size : "+str(np.size(STATE) )
  for idSrc in IDs[FR].keys() :
    for idTrgt in IDs[EN].keys() :
      score = len([h for h in HAPAX.keys() if (IDs[FR][idSrc] in HAPAX[h][FR]) and (IDs[EN][idTrgt] in HAPAX[h][EN]) ])
      #print score
      STATE[idSrc, idTrgt] = float(score)
  filename = "state_"+str(nbIter)
  print filename
  #print STATE
  np.save(filename, STATE)

def reduction() :
  #iterate
  convergence = False
  MAX_ITER = 20
  stop = False
  while not (convergence or nbIter == MAX_ITER or stop) :
    state1 = np.copy(STATE)
    S = np.sum(STATE, axis=0)
    #print S
    S[(S < eps)] = 1
    #print S
    X = np.matlib.repmat(S, nr, 1)
    #print X
    PROB["0"] = STATE / X
    #print PROB["0"]

    #print STATE
    S = np.sum(STATE, axis=1)
    #print S
    S[(S < eps)] = 1
    #print S
    X = np.transpose(np.matlib.repmat(S, nc, 1))
    #print X
    PROB["1"] = STATE / X
    #print PROB["1"]

    #print "\n\n"
    STATE = PROB["0"]*PROB["1"]
    #print P
    #print STATE

    max_col = np.amax(STATE, axis=1)
    print "max = "+str(np.amax(max_col))
    #thres = min(np.amax(max_col)-alpha*eps, thres)
    print "threshold "+str(nbIter)+" = "+str(thres)
    ##print max_col
    for i in range(nc) :
      if max_col[i] > thres :
        #print i
        row = STATE[i,:]
        ind = np.nonzero((row < max_col[i]) & (row > 0.0))
        if len(ind) > 0 :
          print max_col[i]
          STATE[i,ind] = 0.0
        #print STATE
    #print STATE

##    #check invariant : no row is nul

    print "================"
    nbIter = nbIter + 1
    if ( (nbIter % 10) == 0 ) :
      print "Saving state : nbIter = "+str(nbIter)
      filename = "state_"+str(nbIter)
      print filename
      print STATE
      np.save(filename, STATE)
    if ( abs((state1 - STATE)) < np.matlib.repmat(eps, nr, nc) ).all() :
      #print STATE
      convergence = True
      print "Convergence. nbIter="+str(nbIter)
      print float(1)/400000
      print eps
  print "Saving final result..."
  filename = "state_final"
  print filename
  np.save(filename, STATE)

  print "nbIter = "+str(nbIter)+"\t Convergence = "+str(convergence)+"\t stop = "+str(stop)
