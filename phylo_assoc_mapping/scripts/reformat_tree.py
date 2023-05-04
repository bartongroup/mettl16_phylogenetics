import sys
import ete3

# combine the two trees outputted by Orthofinder into a single tree with NHX confidence tags

if __name__ == '__main__':
    tree_named = ete3.Tree(sys.argv[1], format=1)
    tree_support = ete3.Tree(sys.argv[2])
    tree_out = sys.argv[3]
    for node_named, node_support in zip(tree_named.traverse(), tree_support.traverse()):
        assert node_named.dist == node_support.dist, 'Trees do not match'
        if node_named.is_leaf():
            assert node_support.is_leaf(), 'Trees do not match'
            assert node_named.name == node_support.name, 'Trees do not match'
        node_named.add_feature('conf', node_support.support)
    with open(tree_out, 'w') as f:
        f.write(tree_named.write(features=['conf'], format=1, format_root_node=True))