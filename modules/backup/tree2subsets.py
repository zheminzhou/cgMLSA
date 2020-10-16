import sys, click, numpy as np
from ete3 import Tree
from _collections import defaultdict


@click.command()
@click.argument('tree')
@click.argument('group_num', type=int)
def main(tree, group_num) :
    tre = Tree(tree, format=1)
    for nodeId, node in enumerate(tre.traverse('postorder')) :
        if node.is_leaf() :
            node.descendents = {node.name:[0, 0]}
        else :
            node.name = 'N_{0}'.format(nodeId)
            node.descendents = {}
            for c in node.children :
                for desc, (nbr, lbr) in c.descendents.items() :
                    node.descendents[desc] = [nbr+1, lbr+c.dist]

    for nodeId, node in enumerate(tre.traverse('postorder')) :
        if not node.is_leaf() :
            node.descendents2 = {}
            for c in node.children:
                try :
                    for desc, (nbr, lbr) in c.descendents2.items():
                        node.descendents2[desc] = [nbr + 1, lbr + c.dist]
                except :
                    for desc, (nbr, lbr) in c.descendents.items():
                        node.descendents2[desc] = [nbr + 1, lbr + c.dist]

            if len(node.descendents2) >= group_num or not node.up :
                inGroup = [ desc for desc, info in node.descendents2.items()]
                sys.stdout.write('{0}\t{1}'.format(nodeId, ','.join(inGroup)))
                if node.up :
                    outgroup = min([[desc, info ] for desc, info in node.up.descendents.items() \
                                       if desc not in node.descendents], key=lambda x:x[1])
                    sys.stdout.write('\t{0}\n'.format(outgroup[0]))
                else :
                    sys.stdout.write('\t\n'.format(outgroup[0]))

                d = min([[desc, info] for desc, info in node.descendents2.items()], \
                                         key=lambda x:x[1])
                node.descendents2 = {d[0]:d[1]}

if __name__ == '__main__' :
    main()