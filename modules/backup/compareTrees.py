import os, sys, click, numpy as np, pandas as pd, re, subprocess
from ete3 import Tree
from uberBlast import uberBlast
from configure import readFasta, uopen, rc, externals

def delNode(tre, node):
    dist = node.dist
    parent = node.up
    if parent is None :
        if tre != node :
            return tre
        parent = node.children[0]
        dist = parent.dist
        node.remove_child(parent)
        parent.up = None
        tre = parent
    else :
        parent.remove_child(node)
        
    for child in node.children :
        parent.add_child(child)
        child.up = parent
        child.dist += dist
    return tre

def delTip(tre, tip):
    parent = tip.up
    if parent is None :
        return delNode(tre, tip)
    parent.remove_child(tip)
    tip.up = None
    if len(parent.children) < 2 :
        tre = delNode(tre, parent)
    return tre

def getNumber(s, default=0) :
    try :
        return float(s)
    except :
        return default
    

@click.command()
@click.argument('tree1', nargs=1)
@click.argument('tree2', nargs=1)
def main(tree1, tree2) :
    tre1 = Tree(tree1, format=1)
    tre2 = Tree(tree2, format=1)
    shared_tips = set(tre1.get_leaf_names()) & set(tre2.get_leaf_names())
    for tip in tre1.get_leaves() :
        if tip.name not in shared_tips :
            tre1 = delTip(tre1, tip)
    for tip in tre2.get_leaves() :
        if tip.name not in shared_tips :
            tre2 = delTip(tre2, tip)
    with uopen('tre1', label='x') as t1, uopen('tre2', label='x') as t2 :
        t1.close()
        t2.close()
        tre1.write(outfile=t1.name, format=1)
        tre2.write(outfile=t2.name, format=1)
        o1 = subprocess.Popen('{tqDist} -v {0} {1}'.format(t1.name, t1.name, **externals).split(), stdout=subprocess.PIPE).communicate()
        o2 = subprocess.Popen('{tqDist} -v {0} {1}'.format(t2.name, t2.name, **externals).split(), stdout=subprocess.PIPE).communicate()
        o3 = subprocess.Popen('{tqDist} -v {0} {1}'.format(t1.name, t2.name, **externals).split(), stdout=subprocess.PIPE).communicate()
        o1 = [float(x) for x in o1[0].split(b'\t')]
        o2 = [float(x) for x in o2[0].split(b'\t')]
        o3 = [float(x) for x in o3[0].split(b'\t')]
        varied = (o1[6]+o2[6]-2*o3[6])
        dist = o3[2]-varied
        shared = o3[4]
        summed = varied+dist+shared
        print('{0}\t{1}\t{2}\t{3}\t{4}'.format(tree1, tree2, shared/summed, dist/summed, varied/summed))
    
if __name__ == '__main__' :
    main()