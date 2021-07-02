package jneat.utils;

import jneat.Gene;
import jneat.Genome;
import jneat.NNode;

public class GenomeVerifier {
    private final Genome genome;

    public GenomeVerifier(Genome genome) {
        this.genome = genome;
    }

    public boolean verify() {
        if (genome.getGenes().size() == 0 || genome.getNodes().size() == 0 || genome.getTraits().size() == 0) {
            return false;
        }

        // control if nodes in gene are defined and are the same nodes il nodes list
        for (Gene _gene : genome.getGenes()) {
            NNode inode = _gene.getLnk().getIn_node();
            NNode onode = _gene.getLnk().getOut_node();

            if (inode == null) {
                System.out.println(" *ERROR* inode = null in genome #" + genome.genome_id);
                return false;
            }

            if (onode == null) {
                System.out.println(" *ERROR* onode = null in genome #" + genome.genome_id);
                return false;
            }

            if (!genome.getNodes().contains(inode)) {
                System.out.println("Missing inode:  node defined in gene not found in Vector nodes of genome #" + genome.genome_id);
                System.out.print("\n the inode is=" + inode.getNode_id());
                return false;
            }

            if (!genome.getNodes().contains(onode)) {
                System.out.println("Missing onode:  node defined in gene not found in Vector nodes of genome #" + genome.genome_id);
                System.out.print("\n the onode is=" + onode.getNode_id());
                return false;
            }
        }

        // verify if list nodes is ordered
        int last_id = 0;
        for (NNode _node : genome.getNodes()) {
            if (_node.getNode_id() < last_id) {
                System.out.println("ALERT: NODES OUT OF ORDER : ");
                System.out.println(" last node_id is= " + last_id + " , current node_id=" + _node.getNode_id());
                return false;
            }
            last_id = _node.getNode_id();
        }

        // control in genes are gene duplicate for contents
        for (Gene _gene : genome.getGenes()) {
            int i1 = _gene.getLnk().getIn_node().getNode_id();
            int o1 = _gene.getLnk().getOut_node().getNode_id();
            boolean r1 = _gene.getLnk().getIs_recurrent();

            for (Gene _gene1 : genome.getGenes()) {
                if (_gene1.getLnk().getIn_node().getNode_id() == i1
                        && _gene1.getLnk().getOut_node().getNode_id() == o1
                        && _gene1.getLnk().getIs_recurrent() == r1) {
                    System.out.print(" \n  ALERT: DUPLICATE GENES :");
                    System.out.print("  inp_node=" + i1 + " out_node=" + o1);
                    System.out.print("  in GENOME id -->" + genome.genome_id);
                    System.out.print("  gene1 is : ");
                    _gene.op_view();
                    System.out.print("  gene2 is : ");
                    _gene1.op_view();

                    return false;
                }
            }
        }

        if (genome.getNodes().size() >= 500) {
            boolean disab = false;
            for (Gene _gene : genome.getGenes()) {
                if (!_gene.getEnable() && disab) {
                    System.out.print("\n ALERT: 2 DISABLES IN A ROW: " + _gene.getLnk().getIn_node().getNode_id());
                    System.out.print(" inp node=" + _gene.getLnk().getIn_node().getNode_id());
                    System.out.print(" out node=" + _gene.getLnk().getOut_node().getNode_id());
                    System.out.print(" for GENOME " + genome.genome_id);
                    System.out.print("\n Gene is :");
                    _gene.op_view();
                }

                disab = !_gene.getEnable();
            }
        }

        return true;
    }
}
