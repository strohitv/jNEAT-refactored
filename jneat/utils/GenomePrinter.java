package jneat.utils;

import jNeatCommon.NeatConstant;
import jneat.Gene;
import jneat.Genome;
import jneat.NNode;
import jneat.Trait;

public class GenomePrinter {
    private final Genome genome;

    public GenomePrinter(Genome genome) {
        this.genome = genome;
    }

    public void printGenome() {
        System.out.print("\n GENOME START   id=" + genome.genome_id);
        System.out.print("\n  nodes are :" + genome.getNodes().size());
        System.out.print("\n  genes are :" + genome.getGenes().size());
        System.out.print("\n  trait are :" + genome.getTraits().size());

        for (NNode _node : genome.getNodes()) {
            if (_node.getGen_node_label() == NeatConstant.INPUT) {
                System.out.print("\n Input ");
            }

            if (_node.getGen_node_label() == NeatConstant.OUTPUT) {
                System.out.print("\n Output");
            }

            if (_node.getGen_node_label() == NeatConstant.HIDDEN) {
                System.out.print("\n Hidden");
            }

            if (_node.getGen_node_label() == NeatConstant.BIAS) {
                System.out.print("\n Bias  ");
            }

            _node.op_view();
        }

        for (Gene _gene : genome.getGenes()) {
            _gene.op_view();
        }

        System.out.print("\n");
        System.out.print(" Traits:\n");

        for (Trait _trait : genome.getTraits()) {
            _trait.op_view();
        }

        System.out.print("\n");
        System.out.print(" GENOME END");
    }
}
