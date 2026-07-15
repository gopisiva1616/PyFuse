#from pycallgraph import PyCallGraph
#from pycallgraph.output import GraphvizOutput

from pyfuse.cli import main

if __name__ == "__main__":
    #graphviz = GraphvizOutput()
    #graphviz.output_file = 'code_map.png'
    #with PyCallGraph(output=graphviz):
    raise SystemExit(main())