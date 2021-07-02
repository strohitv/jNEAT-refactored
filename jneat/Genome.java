package jneat;

import jNeatCommon.IOseq;
import jNeatCommon.NeatConstant;
import jNeatCommon.NeatRoutine;
import jneat.utils.CompabilityCounter;

import java.text.DecimalFormat;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.stream.Collectors;

public class Genome extends Neat {
    // Is a reference from this genotype to phenotype
    Network phenotype;

    // Numeric identification for this genotype
    public int genome_id;

    // Each Gene in (3) has a marker telling when it arose historically;
    // Thus, these Genes can be used to speciate the population, and
    // the list of Genes provide an evolutionary history of innovation and link-building
    Vector<Gene> genes;

    // parameter conglomerations :Reserved parameter space for future use
    Vector<Trait> traits;

    // Is a collection of NNode mapped in a Vector;
    Vector<NNode> nodes;

    // note are two String for store statistics information
    // when genomes are readed (if exist : null otherwise);
    public String notes;

    public int getGenome_id() {
        return genome_id;
    }

    public Vector<Gene> getGenes() {
        return genes;
    }

    public Vector<Trait> getTraits() {
        return traits;
    }

    public Vector<NNode> getNodes() {
        return nodes;
    }

    public void setNodes(Vector<NNode> nodes) {
        this.nodes = nodes;
    }

    public Network getPhenotype() {
        return phenotype;
    }

    public Genome duplicate(int new_id) {
        Vector<Trait> traits_dup = new Vector<>(traits.size(), 0);
        Vector<NNode> nodes_dup = new Vector<>(nodes.size(), 0);
        Vector<Gene> genes_dup = new Vector<>(genes.size(), 0);

        // duplicate trait
        traits_dup.addAll(traits.stream().map(Trait::new).collect(Collectors.toList()));

        // duplicate NNodes
        for (NNode _node : nodes) {
            NNode newnode = new NNode(
                    _node,
                    traits_dup.stream().filter(t -> _node.nodetrait != null && t.trait_id == _node.nodetrait.trait_id).findFirst().orElse(null)
            );
            _node.dup = newnode;
            nodes_dup.add(newnode);
        }

        // duplicate Genes
        for (Gene _gene : genes) {
            // creation of new gene with a pointer to new node
            genes_dup.add(new Gene(
                    _gene,
                    traits_dup.stream().filter(t -> _gene.lnk.linktrait != null && t.trait_id == _gene.lnk.linktrait.trait_id).findFirst().orElse(null),
                    _gene.lnk.in_node.dup,
                    _gene.lnk.out_node.dup)
            );
        }

        // okay all nodes created, the new genome can be generate
        return new Genome(new_id, traits_dup, nodes_dup, genes_dup);
    }

    public Genome(int id, Vector<Trait> traits, Vector<NNode> nodes, Vector<Gene> genes) {
        genome_id = id;
        this.traits = traits;
        this.nodes = nodes;
        this.genes = genes;
        notes = null;
        phenotype = null;
    }

    public void mutate_link_weight(double power, double rate, int mutation_type) {
        //The power of mutation will rise farther into the genome
        //on the theory that the older genes are more fit since
        //they have stood the test of time
        boolean severe = NeatRoutine.randfloat() > 0.5; //Once in a while really shake things up // for 50% of Prob. severe is true
        double num = 0.0; //counts gene placement

        for (Gene _gene : genes) {
            double gausspoint = 0.3;
            double coldgausspoint = 0.1;

            if (!severe) { // with other 50%.....
                if ((genes.size() >= 10.0) && (num > genes.size() * 0.8)) { // 0.8 => Signifies the last part of the genome
                    gausspoint = 0.5;
                    coldgausspoint = 0.3;
                } else if (NeatRoutine.randfloat() > 0.5) {
                    gausspoint = 1.0 - rate;
                    coldgausspoint = 1.0 - rate - 0.1;
                } else {
                    gausspoint = 1.0 - rate;
                    coldgausspoint = 1.0 - rate;
                }
            }

            // choose a number from ]-1,+1[
            double randnum = NeatRoutine.randposneg() * NeatRoutine.randfloat() * power;

            if (mutation_type == NeatConstant.GAUSSIAN) {
                double randchoice = NeatRoutine.randfloat(); // a number from ]0,1[  //Decide what kind of mutation to do on a gene
                if (randchoice > gausspoint) {
                    _gene.lnk.weight += randnum;
                } else if (randchoice > coldgausspoint) {
                    _gene.lnk.weight = randnum;
                }
            } else if (mutation_type == NeatConstant.COLDGAUSSIAN) {
                _gene.lnk.weight = randnum;
            }

            // copy to mutation_num, the current weight
            _gene.mutation_num = _gene.lnk.weight;
            num += 1.0;
        }
    }

    public Network genesis(int id) {
        Vector<NNode> inlist = new Vector<>(1, 0);
        Vector<NNode> outlist = new Vector<>(1, 0);
        Vector<NNode> all_list = new Vector<>(nodes.size(), 0);

        for (NNode _node : nodes) {
            // create a copy of gene node for phenotype
            NNode newnode = new NNode(_node.type, _node.node_id);
            newnode.derive_trait(_node.nodetrait);
            newnode.inner_level = 0;
            newnode.gen_node_label = _node.gen_node_label;
            newnode.is_traversed = false;

            if (_node.gen_node_label == NeatConstant.INPUT || _node.gen_node_label == NeatConstant.BIAS) {
                inlist.add(newnode);
            }

            if (_node.gen_node_label == NeatConstant.OUTPUT) {
                outlist.add(newnode);
            }

            all_list.add(newnode);
            _node.analogue = newnode;
        }

        if (genes.size() == 0) {
            System.out.print("\n ALERT : are a network whitout GENES; the result can unpredictable");
        }

        if (outlist.size() == 0) {
            System.out.print("\n ALERT : are a network whitout OUTPUTS; the result can unpredictable");
            this.op_view();
        }

        for (Gene _gene : genes.stream().filter(Gene::getEnable).collect(Collectors.toList())) {
            // Only create the link if the gene is enabled
            NNode inode = _gene.lnk.in_node.analogue;
            NNode onode = _gene.lnk.out_node.analogue;

            // NOTE: This line could be run through a recurrency check if desired (no need to in the current implementation of NEAT)
            Link newlink = new Link(_gene.lnk.weight, inode, onode, _gene.lnk.is_recurrent);
            onode.incoming.add(newlink);
            inode.outgoing.add(newlink);

            // Derive link's parameters from its Trait pointer of linktrait
            _gene.lnk.derive_trait(_gene.lnk.linktrait);
        }

        //  Create the new network
        //  Attach genotype and phenotype together: newnet point to owner genotype (this)
        //  genotype point to owner phenotype (newnet)
        phenotype = new Network(inlist, outlist, all_list, id);
        phenotype.setGenotype(this);
        return phenotype;
    }

    /**
     * This function gives a measure of compatibility between
     * two Genomes by computing a linear combination of 3
     * characterizing variables of their compatibilty.
     * The 3 variables represent PERCENT DISJOINT GENES,
     * PERCENT EXCESS GENES, MUTATIONAL DIFFERENCE WITHIN
     * MATCHING GENES.  So the formula for compatibility
     * is:  disjoint_coeff*pdg+excess_coeff*peg+mutdiff_coeff*mdmg.
     * The 3 coefficients are global system parameters
     */
    public double compatibility(Genome other) {
        CompabilityCounter counter = new CompabilityCounter();

        // Get the length of the longest Genome for percentage computations
        double max_genome_size = Math.max(genes.size(), other.genes.size()); //Size of larger Genome

        // Now move through the Genes of each potential parent until both Genomes end
        for (int j = 0; j < max_genome_size; j++) {
            if (counter.getFirstGenePosition() >= genes.size()) {
                counter.addOtherExcess();
            } else if (counter.getSecondGenePosition() >= other.genes.size()) {
                counter.addOwnExcess();
            } else {
                Gene _gene1 = genes.elementAt(counter.getFirstGenePosition());
                Gene _gene2 = other.genes.elementAt(counter.getSecondGenePosition());

                if (_gene1.innovation_num == _gene2.innovation_num) {
                    counter.addMatching(_gene1.mutation_num, _gene2.mutation_num);
                } else if (_gene1.innovation_num < _gene2.innovation_num) {
                    counter.addOwnDisjoint();
                } else {
                    counter.addOtherDisjoint();
                }
            }
        }

        // Return the compatibility number using compatibility formula
        // Note that mut_diff_total/num_matching gives the AVERAGE difference between mutation_nums for any two matching Genes in the Genome.
        // Look at disjointedness and excess in the absolute (ignoring size)
        return (Neat.p_disjoint_coeff * (counter.getDisjointCount() / 1.0) + Neat.p_excess_coeff * (counter.getExcessCount() / 1.0)
                + Neat.p_mutdiff_coeff * (counter.getTotalWeigthDifference() / counter.getMatchingCount()));
    }

    public double get_last_gene_innovnum() {
        return genes.lastElement().innovation_num + 1;
    }

    public int get_last_node_id() {
        return nodes.lastElement().node_id + 1;
    }

    public void op_view() {
        System.out.print("\n GENOME START   id=" + genome_id);
        System.out.print("\n  genes are :" + genes.size());
        System.out.print("\n  nodes are :" + nodes.size());
        System.out.print("\n  trait are :" + traits.size());

        for (NNode _node : nodes) {
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

        for (Gene _gene : genes) {
            _gene.op_view();
        }

        System.out.print("\n");
        System.out.print(" Traits:\n");

        for (Trait _trait : traits) {
            _trait.op_view();
        }

        System.out.print("\n");
        System.out.print(" GENOME END");
    }

    public boolean verify() {
        if (genes.size() == 0 || nodes.size() == 0 || traits.size() == 0) {
            return false;
        }

        // control if nodes in gene are defined and are the same nodes il nodes list
        for (Gene _gene : genes) {
            NNode inode = _gene.lnk.in_node;
            NNode onode = _gene.lnk.out_node;

            if (inode == null) {
                System.out.println(" *ERROR* inode = null in genome #" + genome_id);
                return false;
            }

            if (onode == null) {
                System.out.println(" *ERROR* onode = null in genome #" + genome_id);
                return false;
            }

            if (!nodes.contains(inode)) {
                System.out.println("Missing inode:  node defined in gene not found in Vector nodes of genome #" + genome_id);
                System.out.print("\n the inode is=" + inode.node_id);
                return false;
            }

            if (!nodes.contains(onode)) {
                System.out.println("Missing onode:  node defined in gene not found in Vector nodes of genome #" + genome_id);
                System.out.print("\n the onode is=" + onode.node_id);
                return false;
            }
        }

        // verify if list nodes is ordered
        int last_id = 0;
        for (NNode _node : nodes) {
            if (_node.node_id < last_id) {
                System.out.println("ALERT: NODES OUT OF ORDER : ");
                System.out.println(" last node_id is= " + last_id + " , current node_id=" + _node.node_id);
                return false;
            }
            last_id = _node.node_id;
        }

        // control in genes are gene duplicate for contents
        for (Gene _gene : genes) {
            int i1 = _gene.lnk.in_node.node_id;
            int o1 = _gene.lnk.out_node.node_id;
            boolean r1 = _gene.lnk.is_recurrent;

            for (Gene _gene1 : genes) {
                if (_gene1.lnk.in_node.node_id == i1
                        && _gene1.lnk.out_node.node_id == o1
                        && _gene1.lnk.is_recurrent == r1) {
                    System.out.print(" \n  ALERT: DUPLICATE GENES :");
                    System.out.print("  inp_node=" + i1 + " out_node=" + o1);
                    System.out.print("  in GENOME id -->" + genome_id);
                    System.out.print("  gene1 is : ");
                    _gene.op_view();
                    System.out.print("  gene2 is : ");
                    _gene1.op_view();

                    return false;
                }
            }
        }

        if (nodes.size() >= 500) {
            boolean disab = false;
            for (Gene _gene : genes) {
                if (!_gene.enable && disab) {
                    System.out.print("\n ALERT: 2 DISABLES IN A ROW: " + _gene.lnk.in_node.node_id);
                    System.out.print(" inp node=" + _gene.lnk.in_node.node_id);
                    System.out.print(" out node=" + _gene.lnk.out_node.node_id);
                    System.out.print(" for GENOME " + genome_id);
                    System.out.print("\n Gene is :");
                    _gene.op_view();
                }

                disab = !_gene.enable;
            }
        }

        return true;
    }

    public void print_to_filename(String xNameFile) {
        // write to file genome in native format (for re-read)
        IOseq xFile;

        xFile = new IOseq(xNameFile);
        xFile.IOseqOpenW(false);

        try {
            print_to_file(xFile);
        } catch (Throwable e) {
            e.printStackTrace();
        }

        xFile.IOseqCloseW();
    }

    public Genome mate_multipoint(Genome g, int genomeid, double fitness1, double fitness2) {
        boolean disable = false; //Set to true if we want to disabled a chosen gene
        NNode curnode = null;
        Gene chosengene = null;

        //Tells if the first genome (this one) has better fitness or not
        boolean skip;

        //First, average the Traits from the 2 parents to form the baby's Traits
        //It is assumed that trait vectors are the same length
        //In the future, may decide on a different method for
        //trait mating (corrispondenza)
        Vector<Trait> newtraits = new Vector<>(traits.size(), 0);

        for (int j = 0; j < traits.size(); j++) {
            newtraits.add(new Trait(traits.elementAt(j), g.traits.elementAt(j)));
        }

        //Figure out which genome is better
        //The worse genome should not be allowed to add extra structural baggage
        //If they are the same, use the smaller one's disjoint and excess genes only
        boolean p1better = false;

        int size1 = genes.size();
        int size2 = g.genes.size();

        if (fitness1 > fitness2) {
            p1better = true;
        } else if (fitness1 == fitness2 && size1 < size2) {
            p1better = true;
        }

        int len_genome = Math.max(size1, size2);
        int len_nodes = nodes.size();

        Vector<Gene> newgenes = new Vector<>(len_genome, 0);
        Vector<NNode> newnodes = new Vector<>(len_nodes, 0);

        int j1 = 0;
        int j2 = 0;
        while (j1 < size1 || j2 < size2) {
            //  chosen of 'just' gene
            skip = false; //Default to not skipping a chosen gene
            if (j1 >= size1) {
                chosengene = g.genes.elementAt(j2);
                j2++;
                if (p1better) {
                    skip = true; //Skip excess from the worse genome
                }
            } else if (j2 >= size2) {
                chosengene = genes.elementAt(j1);
                j1++;
                if (!p1better) {
                    skip = true; //Skip excess from the worse genome
                }
            } else {
                Gene _p1gene = genes.elementAt(j1);
                Gene _p2gene = g.genes.elementAt(j2);

                double p1innov = _p1gene.innovation_num;
                double p2innov = _p2gene.innovation_num;
                if (p1innov == p2innov) {
                    if (NeatRoutine.randfloat() < 0.5) {
                        chosengene = _p1gene;
                    } else {
                        chosengene = _p2gene;
                    }

                    //If one is disabled, the corresponding gene in the offspring
                    //will likely be disabled
                    disable = false;
                    if (!_p1gene.enable || !_p2gene.enable) {
                        if (NeatRoutine.randfloat() < 0.75) {
                            disable = true;
                        }
                    }

                    j1++;
                    j2++;
                } else if (p1innov < p2innov) {
                    chosengene = _p1gene;
                    j1++;
                    if (!p1better) {
                        skip = true;
                    }
                } else if (p2innov < p1innov) {
                    chosengene = _p2gene;
                    j2++;
                    if (p1better) {
                        skip = true;
                    }
                }
            }

            assert chosengene != null;

            //Check to see if the chosengene conflicts with an already chosen gene
            //i.e. do they represent the same link
            for (Gene gene : newgenes) {
                if (gene.lnk.in_node.node_id == chosengene.lnk.in_node.node_id
                        && gene.lnk.out_node.node_id == chosengene.lnk.out_node.node_id
                        && gene.lnk.is_recurrent == chosengene.lnk.is_recurrent) {
                    skip = true;
                    break;
                }

                if (gene.lnk.in_node.node_id == chosengene.lnk.out_node.node_id
                        && gene.lnk.out_node.node_id == chosengene.lnk.in_node.node_id
                        && !gene.lnk.is_recurrent
                        && !chosengene.lnk.is_recurrent) {
                    skip = true;
                    break;
                }
            }

            if (!skip) {
                NNode new_inode;
                NNode new_onode;

                //Now add the chosengene to the baby
                //First, get the trait pointer
                int first_traitnum = traits.firstElement().trait_id;

                int traitnum = first_traitnum;
                if (chosengene.lnk.linktrait != null) {
                    traitnum = chosengene.lnk.linktrait.trait_id - first_traitnum;
                }

                //Next check for the nodes, add them if not in the baby Genome already

                NNode inode = chosengene.lnk.in_node;
                NNode onode = chosengene.lnk.out_node;

                //Check for inode in the newnodes list
                //Check for inode, onode in the newnodes list
                boolean found;
                if (inode.node_id < onode.node_id) {
                    // search the inode
                    found = false;
                    for (int ix = 0; ix < newnodes.size(); ix++) {
                        curnode = newnodes.elementAt(ix);
                        if (curnode.node_id == inode.node_id) {
                            found = true;
                            break;
                        }
                    }

                    // if exist , point to exitsting version
                    if (found) {
                        new_inode = curnode;
                    } else { // else create the inode
                        int nodetraitnum = 0;
                        if (inode.nodetrait != null) {
                            nodetraitnum = inode.nodetrait.trait_id - first_traitnum;
                        }

                        new_inode = new NNode(inode, newtraits.elementAt(nodetraitnum));
                        //insert in newnodes list
                        node_insert(newnodes, new_inode);
                    }

                    // search the onode
                    found = false;
                    for (int ix = 0; ix < newnodes.size(); ix++) {
                        curnode = newnodes.elementAt(ix);
                        if (curnode.node_id == onode.node_id) {
                            found = true;
                            break;
                        }
                    }

                    // if exist , point to exitsting version
                    if (found) {
                        new_onode = curnode;
                    } else { // else create the onode
                        int nodetraitnum = 0;
                        if (onode.nodetrait != null) {
                            nodetraitnum = onode.nodetrait.trait_id - first_traitnum;
                        }

                        new_onode = new NNode(onode, newtraits.elementAt(nodetraitnum));
                        //insert in newnodes list
                        node_insert(newnodes, new_onode);
                    }
                } else { // end block : inode.node_id < onode.node_id
                    // search the onode
                    found = false;
                    for (int ix = 0; ix < newnodes.size(); ix++) {
                        curnode = newnodes.elementAt(ix);
                        if (curnode.node_id == onode.node_id) {
                            found = true;
                            break;
                        }
                    }

                    // if exist , point to exitsting version
                    if (found) {
                        new_onode = curnode;
                    } else { // else create the onode
                        int nodetraitnum = 0;
                        if (onode.nodetrait != null) {
                            nodetraitnum = onode.nodetrait.trait_id - first_traitnum;
                        }

                        new_onode = new NNode(onode, newtraits.elementAt(nodetraitnum));
                        //insert in newnodes list
                        node_insert(newnodes, new_onode);
                    }

                    // search the inode
                    found = false;
                    for (int ix = 0; ix < newnodes.size(); ix++) {
                        curnode = newnodes.elementAt(ix);
                        if (curnode.node_id == inode.node_id) {
                            found = true;
                            break;
                        }
                    }

                    // if exist , point to exitsting version
                    if (found) {
                        new_inode = curnode;
                    } else { // else create the inode
                        int nodetraitnum = 0;
                        if (inode.nodetrait != null) {
                            nodetraitnum = inode.nodetrait.trait_id - first_traitnum;
                        }

                        new_inode = new NNode(inode, newtraits.elementAt(nodetraitnum));
                        //insert in newnodes list
                        node_insert(newnodes, new_inode);
                    }
                }

                //Add the Gene
                Gene newgene = new Gene(chosengene, newtraits.elementAt(traitnum), new_inode, new_onode);
                if (disable) {
                    newgene.enable = false;
                    disable = false;
                }
                newgenes.add(newgene);
            }
        } // end block genome (while)

        Genome new_genome = new Genome(genomeid, newtraits, newnodes, newgenes);

        //	boolean h = new_genome.verify();
        boolean found = false;
        for (int ix = 0; ix < newnodes.size(); ix++) {
            curnode = newnodes.elementAt(ix);
            if (curnode.gen_node_label == NeatConstant.OUTPUT) {
                found = true;
                break;
            }
        }

        if (!found) {
            System.out.print("\n *--------------- not found output node ----------------------------");
            System.out.print("\n * during mate_multipoint : please control the following's *********");
            System.out.print("\n * control block : ");
            System.out.print("\n Genome A= ");
            this.op_view();
            System.out.print("\n Genome B= ");
            g.op_view();
            System.out.print("\n Result = ");
            new_genome.op_view();
            System.exit(0);
        }

        return new_genome;
    }

    public Genome mate_multipoint_avg(Genome g, int genomeid, double fitness1, double fitness2) {
        boolean disable = false; //Set to true if we want to disabled a chosen gene
        int traitnum;

        //Set up the avgene
        Gene avgene = new Gene(null, 0.0, null, null, false, 0.0, 0.0);

        //First, average the Traits from the 2 parents to form the baby's Traits
        //It is assumed that trait vectors are the same length
        //In the future, may decide on a different method for
        //trait mating (corrispondenza)
        int len_traits = traits.size();

        Vector<Trait> newtraits = new Vector<>(len_traits, 0);
        for (int j = 0; j < len_traits; j++) {
            newtraits.add(new Trait(traits.elementAt(j), g.traits.elementAt(j)));
        }

        //Figure out which genome is better
        //The worse genome should not be allowed to add extra structural baggage
        //If they are the same, use the smaller one's disjoint and excess genes only
        boolean p1better = false;

        if (fitness1 > fitness2 || (fitness1 == fitness2 && genes.size() < g.genes.size())) {
            p1better = true;
        }

        int len_genome = Math.max(genes.size(), g.genes.size());
        int len_nodes = nodes.size();

        Vector<Gene> newgenes = new Vector<>(len_genome, 0);
        Vector<NNode> newnodes = new Vector<>(len_nodes, 0);

        Gene chosengene = null;
        NNode curnode = null;

        boolean skip;

        int j1 = 0;
        int j2 = 0;
        while (j1 < genes.size() || j2 < g.genes.size()) {
            //  chosen of 'just' gene
            avgene.enable = true; //Default to enabled
            skip = false; //Default to not skipping a chosen gene

            if (j1 >= genes.size()) {
                chosengene = g.genes.elementAt(j2);
                j2++;
                if (p1better) {
                    skip = true; //Skip excess from the worse genome
                }
            } else if (j2 >= g.genes.size()) {
                chosengene = genes.elementAt(j1);
                j1++;
                if (!p1better) {
                    skip = true; //Skip excess from the worse genome
                }
            } else {
                Gene _p1gene = genes.elementAt(j1);
                Gene _p2gene = g.genes.elementAt(j2);

                if (_p1gene.innovation_num == _p2gene.innovation_num) {
                    if (NeatRoutine.randfloat() > 0.5) {
                        avgene.lnk.linktrait = _p1gene.lnk.linktrait;
                    } else {
                        avgene.lnk.linktrait = _p2gene.lnk.linktrait;
                    }

                    //WEIGHTS AVERAGED HERE
                    avgene.lnk.weight = (_p1gene.lnk.weight + _p2gene.lnk.weight) / 2.0;

                    if (NeatRoutine.randfloat() > 0.5) {
                        avgene.lnk.in_node = _p1gene.lnk.in_node;
                    } else {
                        avgene.lnk.in_node = _p2gene.lnk.in_node;
                    }

                    if (NeatRoutine.randfloat() > 0.5) {
                        avgene.lnk.out_node = _p1gene.lnk.out_node;
                    } else {
                        avgene.lnk.out_node = _p2gene.lnk.out_node;
                    }

                    if (NeatRoutine.randfloat() > 0.5) {
                        avgene.lnk.is_recurrent = _p1gene.lnk.is_recurrent;
                    } else {
                        avgene.lnk.is_recurrent = _p2gene.lnk.is_recurrent;
                    }

                    avgene.innovation_num = _p1gene.innovation_num;
                    avgene.mutation_num = (_p1gene.mutation_num + _p2gene.mutation_num) / 2.0;

                    //If one is disabled, the corresponding gene in the offspring
                    //will likely be disabled
                    disable = false;
                    if (!_p1gene.enable || !_p2gene.enable) {
                        if (NeatRoutine.randfloat() < 0.75) {
                            disable = true;
                        }
                    }

                    chosengene = avgene;

                    j1++;
                    j2++;
                } else if (_p1gene.innovation_num < _p2gene.innovation_num) {
                    chosengene = _p1gene;
                    j1++;
                    if (!p1better) {
                        skip = true;
                    }
                } else if (_p2gene.innovation_num < _p1gene.innovation_num) {
                    chosengene = _p2gene;
                    j2++;
                    if (p1better) {
                        skip = true;
                    }
                }
            } // end chosen gene

            assert chosengene != null;

            //Check to see if the chosengene conflicts with an already chosen gene
            //i.e. do they represent the same link
            for (Gene gene : newgenes) {
                if (gene.lnk.in_node.node_id == chosengene.lnk.in_node.node_id
                        && gene.lnk.out_node.node_id == chosengene.lnk.out_node.node_id
                        && gene.lnk.is_recurrent == chosengene.lnk.is_recurrent) {
                    skip = true;
                    break;
                }

                if (gene.lnk.in_node.node_id == chosengene.lnk.out_node.node_id
                        && gene.lnk.out_node.node_id == chosengene.lnk.in_node.node_id
                        && !gene.lnk.is_recurrent
                        && !chosengene.lnk.is_recurrent) {
                    skip = true;
                    break;
                }
            }

            if (!skip) {
                //Now add the chosengene to the baby
                //First, get the trait pointer
                int first_traitnum = traits.firstElement().trait_id;

                if (chosengene.lnk.linktrait == null) {
                    traitnum = first_traitnum;
                } else {
                    traitnum = chosengene.lnk.linktrait.trait_id - first_traitnum;
                }

                //Next check for the nodes, add them if not in the baby Genome already
                NNode inode = chosengene.lnk.in_node;
                NNode onode = chosengene.lnk.out_node;

                //Check for inode in the newnodes list
                boolean found;
                NNode new_inode;
                NNode new_onode;
                if (inode.node_id < onode.node_id) {
                    // search the inode
                    found = false;
                    for (int ix = 0; ix < newnodes.size(); ix++) {
                        curnode = newnodes.elementAt(ix);
                        if (curnode.node_id == inode.node_id) {
                            found = true;
                            break;
                        }
                    }

                    // if exist , point to exitsting version
                    if (found) {
                        new_inode = curnode;
                    } else { // else create the inode
                        int nodetraitnum = 0;
                        if (inode.nodetrait != null) {
                            nodetraitnum = inode.nodetrait.trait_id - first_traitnum;
                        }

                        new_inode = new NNode(inode, newtraits.elementAt(nodetraitnum));

                        //insert in newnodes list
                        node_insert(newnodes, new_inode);
                    }

                    // search the onode
                    found = false;
                    for (int ix = 0; ix < newnodes.size(); ix++) {
                        curnode = newnodes.elementAt(ix);
                        if (curnode.node_id == onode.node_id) {
                            found = true;
                            break;
                        }
                    }

                    // if exist , point to exitsting version
                    if (found) {
                        new_onode = curnode;
                    } else { // else create the onode
                        int nodetraitnum = 0;
                        if (onode.nodetrait != null) {
                            nodetraitnum = onode.nodetrait.trait_id - first_traitnum;
                        }

                        new_onode = new NNode(onode, newtraits.elementAt(nodetraitnum));

                        //insert in newnodes list
                        node_insert(newnodes, new_onode);
                    }
                } else {// end block : inode.node_id < onode.node_id
                    // search the onode
                    found = false;
                    for (int ix = 0; ix < newnodes.size(); ix++) {
                        curnode = newnodes.elementAt(ix);
                        if (curnode.node_id == onode.node_id) {
                            found = true;
                            break;
                        }
                    }

                    // if exist , point to exitsting version
                    if (found) {
                        new_onode = curnode;
                    } else { // else create the onode
                        int nodetraitnum = 0;
                        if (onode.nodetrait != null) {
                            nodetraitnum = onode.nodetrait.trait_id - first_traitnum;
                        }

                        new_onode = new NNode(onode, newtraits.elementAt(nodetraitnum));

                        //insert in newnodes list
                        node_insert(newnodes, new_onode);
                    }

                    // search the inode
                    found = false;
                    for (int ix = 0; ix < newnodes.size(); ix++) {
                        curnode = newnodes.elementAt(ix);
                        if (curnode.node_id == inode.node_id) {
                            found = true;
                            break;
                        }
                    }

                    // if exist , point to exitsting version
                    if (found) {
                        new_inode = curnode;
                    } else { // else create the inode
                        int nodetraitnum = 0;
                        if (inode.nodetrait != null) {
                            nodetraitnum = inode.nodetrait.trait_id - first_traitnum;
                        }

                        new_inode = new NNode(inode, newtraits.elementAt(nodetraitnum));

                        //insert in newnodes list
                        node_insert(newnodes, new_inode);
                    }
                }

                //Add the Gene
                Gene newgene = new Gene(chosengene, newtraits.elementAt(traitnum), new_inode, new_onode);

                if (disable) {
                    newgene.enable = false;
                    disable = false;
                }

                newgenes.add(newgene);
            }
        } // end block genome

        return new Genome(genomeid, newtraits, newnodes, newgenes);
    }

    public Genome mate_singlepoint(Genome g, int genomeid) {
        Gene chosengene = null;
        int crosspoint;

        NNode curnode = null;

        int traitnum;

        Vector<Trait> newtraits = new Vector<>(traits.size(), 0);

        for (int j = 0; j < traits.size(); j++) {
            newtraits.add(new Trait(traits.elementAt(j), g.traits.elementAt(j)));
        }

        //Set up the avgene
        Gene avgene = new Gene(null, 0.0, null, null, false, 0.0, 0.0);

        Vector<Gene> newgenes = new Vector<>(genes.size(), 0);
        Vector<NNode> newnodes = new Vector<>(nodes.size(), 0);

        int stopA;
        int stopB;
        Vector<Gene> genomeA;
        Vector<Gene> genomeB;
        if (genes.size() < g.genes.size()) {
            crosspoint = NeatRoutine.randint(0, genes.size() - 1);
            stopA = genes.size();
            stopB = g.genes.size();
            genomeA = genes;
            genomeB = g.genes;
        } else {
            crosspoint = NeatRoutine.randint(0, g.genes.size() - 1);
            stopA = g.genes.size();
            stopB = genes.size();
            genomeA = g.genes;
            genomeB = genes;
        }

        boolean done = false;
        double v1 = 0.0;
        double v2 = 0.0;
        double cellA = 0.0;
        double cellB = 0.0;

        int j1 = 0;
        int j2 = 0;

        // compute what is the hight innovation
        double last_innovB = genomeB.elementAt(stopB - 1).innovation_num;
        double cross_innov = 0;

        Gene geneA = null;
        Gene geneB = null;
        int genecounter = 0; //Ready to count to crosspoint
        while (!done) {
            boolean doneA = false;
            boolean doneB = false;
            boolean skip = false; //Default to not skip a Gene
            avgene.enable = true; //Default to true

            if (j1 < stopA) {
                geneA = genomeA.elementAt(j1);
                v1 = geneA.innovation_num;
                doneA = true;
            }

            if (j2 < stopB) {
                geneB = genomeB.elementAt(j2);
                v2 = geneB.innovation_num;
                doneB = true;
            }

            if (doneA && doneB) {
                //
                if (v1 < v2) {
                    cellA = v1;
                    cellB = 0.0;
                    j1++;
                } else if (v1 == v2) {
                    cellA = v1;
                    cellB = v1;
                    j1++;
                    j2++;
                } else {
                    cellA = 0.0;
                    cellB = v2;
                    j2++;
                }
            } else {
                if (doneA) {
                    cellA = v1;
                    cellB = 0.0;
                    j1++;
                } else if (doneB) {
                    cellA = 0.0;
                    cellB = v2;
                    j2++;
                } else {
                    done = true;
                }
            }

            assert geneA != null;
            assert geneB != null;
            if (!done) {
                // innovA = innovB
                if (cellA == cellB) {
                    if (genecounter < crosspoint) {
                        chosengene = geneA;
                        genecounter++;
                    } else if (genecounter == crosspoint) {
                        if (NeatRoutine.randfloat() > 0.5) {
                            avgene.lnk.linktrait = geneA.lnk.linktrait;
                        } else {
                            avgene.lnk.linktrait = geneB.lnk.linktrait;
                        }

                        //WEIGHTS AVERAGED HERE
                        avgene.lnk.weight = (geneA.lnk.weight + geneB.lnk.weight) / 2.0;

                        if (NeatRoutine.randfloat() > 0.5) {
                            avgene.lnk.in_node = geneA.lnk.in_node;
                        } else {
                            avgene.lnk.in_node = geneB.lnk.in_node;
                        }

                        if (NeatRoutine.randfloat() > 0.5) {
                            avgene.lnk.out_node = geneA.lnk.out_node;
                        } else {
                            avgene.lnk.out_node = geneB.lnk.out_node;
                        }

                        if (NeatRoutine.randfloat() > 0.5) {
                            avgene.lnk.is_recurrent = geneA.lnk.is_recurrent;
                        } else {
                            avgene.lnk.is_recurrent = geneB.lnk.is_recurrent;
                        }

                        avgene.innovation_num = geneA.innovation_num;
                        avgene.mutation_num = (geneA.mutation_num + geneB.mutation_num) / 2.0;

                        //If one is disabled, the corresponding gene in the offspring
                        //will likely be disabled

                        if (!geneA.enable || !geneB.enable) {
                            avgene.enable = false;
                        }

                        chosengene = avgene;
                        genecounter++;
                        cross_innov = cellA;
                    } else if (genecounter > crosspoint) {
                        chosengene = geneB;
                        genecounter++;
                    }
                } else if (cellA != 0 && cellB == 0) {
                    // innovA < innovB
                    if (genecounter < crosspoint) {
                        chosengene = geneA; //make geneA
                        genecounter++;
                    } else if (genecounter == crosspoint) {
                        chosengene = geneA;
                        genecounter++;
                        cross_innov = cellA;
                    } else {
                        if (cross_innov > last_innovB) {
                            chosengene = geneA;
                            genecounter++;
                        } else {
                            skip = true;
                        }
                    }
                } else {
                    if (cellA == 0 && cellB != 0) {
                        // innovA > innovB
                        if (genecounter < crosspoint) {
                            skip = true; //skip geneB
                        } else if (genecounter == crosspoint) {
                            skip = true; //skip an illogic case
                        } else {
                            if (cross_innov > last_innovB) {
                                chosengene = geneA; //make geneA
                                genecounter++;
                            } else {
                                chosengene = geneB; //make geneB : this is a pure case o single crossing
                                genecounter++;
                            }
                        }
                    }
                }

                assert chosengene != null;

                for (Gene gene : newgenes) {
                    if (gene.lnk.in_node.node_id == chosengene.lnk.in_node.node_id
                            && gene.lnk.out_node.node_id == chosengene.lnk.out_node.node_id
                            && gene.lnk.is_recurrent == chosengene.lnk.is_recurrent) {
                        skip = true;
                        break;
                    }

                    if (gene.lnk.in_node.node_id == chosengene.lnk.out_node.node_id
                            && gene.lnk.out_node.node_id == chosengene.lnk.in_node.node_id
                            && !gene.lnk.is_recurrent
                            && !chosengene.lnk.is_recurrent) {
                        skip = true;
                        break;
                    }

                } // and else for control of position in gennomeA/B

                if (!skip) {
                    //Now add the chosengene to the baby
                    //First, get the trait pointer
                    int first_traitnum = traits.firstElement().trait_id;

                    if (chosengene.lnk.linktrait == null)
                        traitnum = first_traitnum;
                    else
                        traitnum = chosengene.lnk.linktrait.trait_id - first_traitnum;

                    //Next check for the nodes, add them if not in the baby Genome already

                    NNode inode = chosengene.lnk.in_node;
                    NNode onode = chosengene.lnk.out_node;

                    //Check for inode, onode in the newnodes list
                    boolean found;
                    NNode new_inode;
                    NNode new_onode;
                    if (inode.node_id < onode.node_id) {
                        // search the inode
                        found = false;
                        for (int ix = 0; ix < newnodes.size(); ix++) {
                            curnode = newnodes.elementAt(ix);
                            if (curnode.node_id == inode.node_id) {
                                found = true;
                                break;
                            }
                        }

                        // if exist , point to exitsting version
                        if (found) {
                            new_inode = curnode;
                        } else { // else create the inode
                            int nodetraitnum = 0;
                            if (inode.nodetrait != null) {
                                nodetraitnum = inode.nodetrait.trait_id - first_traitnum;
                            }

                            new_inode = new NNode(inode, newtraits.elementAt(nodetraitnum));

                            //insert in newnodes list
                            node_insert(newnodes, new_inode);
                        }

                        // search the onode
                        found = false;
                        for (int ix = 0; ix < newnodes.size(); ix++) {
                            curnode = newnodes.elementAt(ix);
                            if (curnode.node_id == onode.node_id) {
                                found = true;
                                break;
                            }
                        }

                        // if exist , point to exitsting version
                        if (found) {
                            new_onode = curnode;
                        } else { // else create the onode
                            int nodetraitnum = 0;
                            if (onode.nodetrait != null) {
                                nodetraitnum = onode.nodetrait.trait_id - first_traitnum;
                            }

                            new_onode = new NNode(onode, newtraits.elementAt(nodetraitnum));

                            //insert in newnodes list
                            node_insert(newnodes, new_onode);
                        }
                    } else { // end block : inode.node_id < onode.node_id
                        // search the onode
                        found = false;
                        for (int ix = 0; ix < newnodes.size(); ix++) {
                            curnode = newnodes.elementAt(ix);
                            if (curnode.node_id == onode.node_id) {
                                found = true;
                                break;
                            }
                        }
                        // if exist , point to exitsting version
                        if (found) {
                            new_onode = curnode;
                        } else { // else create the onode
                            int nodetraitnum = 0;
                            if (onode.nodetrait != null) {
                                nodetraitnum = onode.nodetrait.trait_id - first_traitnum;
                            }

                            new_onode = new NNode(onode, newtraits.elementAt(nodetraitnum));

                            //insert in newnodes list
                            node_insert(newnodes, new_onode);
                        }

                        // search the inode
                        found = false;
                        for (int ix = 0; ix < newnodes.size(); ix++) {
                            curnode = newnodes.elementAt(ix);
                            if (curnode.node_id == inode.node_id) {
                                found = true;
                                break;
                            }
                        }

                        // if exist , point to exitsting version
                        if (found) {
                            new_inode = curnode;
                        } else { // else create the inode
                            int nodetraitnum = 0;
                            if (inode.nodetrait != null) {
                                nodetraitnum = inode.nodetrait.trait_id - first_traitnum;
                            }

                            new_inode = new NNode(inode, newtraits.elementAt(nodetraitnum));

                            //insert in newnodes list
                            node_insert(newnodes, new_inode);
                        }
                    }

                    //Add the Gene
                    newgenes.add(new Gene(chosengene, newtraits.elementAt(traitnum), new_inode, new_onode));
                } // end of block gene creation if !skip
            }
        }

        // search the existence of output node
        // if no dump
        return new Genome(genomeid, newtraits, newnodes, newgenes);
    }

    public void mutate_gene_reenable() {
        Gene _gene;

        for (Gene gene : genes) {
            _gene = gene;
            if (!_gene.enable) {
                _gene.enable = true;
                break;
            }
        }
    }

    /**
     * This chooses a random gene, extracts the link from it,
     * and repoints the link to a random trait
     */
    public void mutate_link_trait(int times) {
        int loop;
        Gene _gene;
        Trait _trait;

        for (loop = 1; loop <= times; loop++) {
            _gene = genes.elementAt(NeatRoutine.randint(0, genes.size() - 1));
            _trait = traits.elementAt(NeatRoutine.randint(0, (traits.size()) - 1));
            _gene.lnk.linktrait = _trait;
        }
    }

    /**
     * This chooses a random node
     * and repoints the node to a random trait
     */
    public void mutate_node_trait(int times) {
        int loop;
        NNode _node;

        for (loop = 1; loop <= times; loop++) {
            _node = nodes.elementAt(NeatRoutine.randint(0, nodes.size() - 1));
            _node.nodetrait = traits.elementAt(NeatRoutine.randint(0, (traits.size()) - 1));
        }
    }

    public void mutate_random_trait() {
        //Retrieve the trait and mutate it
        Trait _trait = traits.elementAt(NeatRoutine.randint(0, (traits.size()) - 1));
        _trait.mutate();
    }

    //
    // Toggle genes from enable on to enable off or
    //   vice versa.  Do it times times.
    public void mutate_toggle_enable(int times) {
        for (int count = 1; count <= times; count++) {
            //find a random gene
            Gene _gene = genes.elementAt(NeatRoutine.randint(0, genes.size() - 1));

            //Toggle the enable on this gene
            if (_gene.enable) {
                //We need to make sure that another gene connects out of the in-node
                //Because if not a section of network will break off and become isolated
                boolean done = false;
                for (int j = 0; j < genes.size(); j++) {
                    Gene _jgene = genes.elementAt(j);
                    if ((_gene.lnk.in_node == _jgene.lnk.in_node)
                            && _jgene.enable
                            && (_jgene.innovation_num != _gene.innovation_num)) {
                        done = true;
                        break;
                    }
                }

                //Disable the gene if it's safe to do so
                if (done) {
                    _gene.enable = false;
                }
            } else {
                _gene.enable = true;
            }
        }
    }

    public void node_insert(Vector<NNode> nlist, NNode n) {
        int position = 0;
        for (int j = 0; j < nlist.size(); j++) {
            if (nlist.elementAt(j).node_id >= n.node_id) {
                position = j;
                break;
            }
        }

        nlist.insertElementAt(n, position);
    }

    public void mutate_add_link(Population pop, int tries) {
        //Find the first non-sensor so that the to-node won't look at sensors as
        //possible destinations
        int first_nonsensor = 0;

        for (NNode node : nodes) {
            if (node.type != NeatConstant.SENSOR) {
                break;
            }
            first_nonsensor++;
        }

        //Decide whether to make this recurrent
        boolean do_recur = NeatRoutine.randfloat() < Neat.p_recur_only_prob;

        boolean found = false;
        NNode thenode1 = null;
        NNode thenode2 = null;

        int trycount = 0;
        while (trycount < tries) {
            // recurrency case .........
            int nodenum1;
            int nodenum2;
            if (do_recur) {
                // at this point :
                // 50% of prob to decide a loop recurrency( node X to node X)
                // 50% a normal recurrency ( node X to node Y)
                if (NeatRoutine.randfloat() > 0.5) {
                    nodenum1 = NeatRoutine.randint(first_nonsensor, nodes.size() - 1);
                    nodenum2 = nodenum1;
                } else {
                    nodenum1 = NeatRoutine.randint(0, nodes.size() - 1);
                    nodenum2 = NeatRoutine.randint(first_nonsensor, nodes.size() - 1);
                }
            } else { // no recurrency case .........
                nodenum1 = NeatRoutine.randint(0, nodes.size() - 1);
                nodenum2 = NeatRoutine.randint(first_nonsensor, nodes.size() - 1);
            }

            // now point to object's nodes
            thenode1 = nodes.elementAt(nodenum1);
            thenode2 = nodes.elementAt(nodenum2);

            // verify if the possible new gene already EXIST
            boolean bypass = false;
            for (int j = 0; j < genes.size(); j++) {
                Gene _gene = genes.elementAt(j);
                if (thenode2.type == NeatConstant.SENSOR) {
                    bypass = true;
                    break;
                }

                if (_gene.lnk.in_node == thenode1
                        && _gene.lnk.out_node == thenode2
                        && _gene.lnk.is_recurrent
                        && do_recur) {
                    bypass = true;
                    break;
                }

                if (_gene.lnk.in_node == thenode1
                        && _gene.lnk.out_node == thenode2
                        && !_gene.lnk.is_recurrent
                        && !do_recur) {
                    bypass = true;
                    break;
                }
            }

            if (!bypass) {
                phenotype.status = 0;
                boolean recurflag = phenotype.has_a_path(thenode1.analogue, thenode2.analogue, 0, nodes.size() * nodes.size());

                if (phenotype.status == 8) {
                    System.out.println("\n  network.mutate_add_link : LOOP DETECTED DURING A RECURRENCY CHECK");
                    return;
                }

                if ((!recurflag && do_recur) || (recurflag && !do_recur)) {
                    trycount++;
                } else {
                    trycount = tries;
                    found = true;
                }
            } // end block bypass

            // if bypass is true, this gene is not good
            // and skip to next cycle
            else {
                trycount++;
            }
        } // end block trycount

        if (found) {
            Gene new_gene = null;
            //Check to see if this innovation already occured in the population
            for (int i = 0; i < pop.innovations.size(); i++) {
                Innovation _innov = pop.innovations.get(i);
                if ((_innov.innovation_type == NeatConstant.NEWLINK)
                        && (_innov.node_in_id == thenode1.node_id)
                        && (_innov.node_out_id == thenode2.node_id)
                        && (_innov.recur_flag == do_recur)) {

                    new_gene = new Gene(traits.elementAt(_innov.new_traitnum), _innov.new_weight, thenode1, thenode2, do_recur, _innov.innovation_num1, 0);
                    break;
                }
            }

            if (new_gene == null) {
                //Choose a random trait
                int traitnum = NeatRoutine.randint(0, traits.size() - 1);

                //Choose the new weight
                double new_weight = NeatRoutine.randposneg() * NeatRoutine.randfloat() * 10.0;

                // read from population current innovation value
                // read curr innovation with postincrement
                double curr_innov = pop.getCurr_innov_num_and_increment();
                //Create the new gene
                new_gene = new Gene(traits.elementAt(traitnum), new_weight, thenode1, thenode2, do_recur, curr_innov, new_weight);
                //Add the innovation
                pop.innovations.add(new Innovation(thenode1.node_id, thenode2.node_id, curr_innov, new_weight, traitnum));
            }

            genes.add(new_gene);
        }
    }

    public void mutate_add_node(Population pop) {
        Gene _gene = null;
        boolean found = false;

        if (genes.size() < 15) {
            boolean step2 = false;
            for (int j = 0; j < genes.size(); j++) {
                _gene = genes.elementAt(j);
                if (_gene.enable && (_gene.lnk.in_node.gen_node_label != NeatConstant.BIAS)) {
                    break;
                }
            }

            for (int j = 0; j < genes.size(); j++) {
                _gene = genes.elementAt(j);
                if ((NeatRoutine.randfloat() >= 0.3) && (_gene.lnk.in_node.gen_node_label != NeatConstant.BIAS)) {
                    step2 = true;
                    break;
                }
            }

            found = step2 && _gene.enable;
        } else {
            for (int trycount = 0; trycount < 20; trycount++) {
                _gene = genes.elementAt(NeatRoutine.randint(0, genes.size() - 1));
                if (_gene.enable && (_gene.lnk.in_node.gen_node_label != NeatConstant.BIAS)) {
                    found = true;
                    break;
                }
            }
        }

        if (!found) {
            return;
        }

        _gene.enable = false;

        Gene newgene1 = null;
        Gene newgene2 = null;
        NNode new_node = null;

        for (int i = 0; i < pop.innovations.size(); i++) {
            Innovation _innov = pop.innovations.get(i);

            if ((_innov.innovation_type == NeatConstant.NEWNODE)
                    && (_innov.node_in_id == _gene.lnk.in_node.node_id)
                    && (_innov.node_out_id == _gene.lnk.out_node.node_id)
                    && (_innov.old_innov_num == _gene.innovation_num)) {
                // Create the new Genes
                // pass this current nodeid to newnode
                new_node = new NNode(NeatConstant.NEURON, _innov.newnode_id, NeatConstant.HIDDEN);
                new_node.nodetrait = traits.firstElement();

                newgene1 = new Gene(_gene.lnk.linktrait, 1.0, _gene.lnk.in_node, new_node, _gene.lnk.is_recurrent, _innov.innovation_num1, 0);
                newgene2 = new Gene(_gene.lnk.linktrait, _gene.lnk.weight, new_node, _gene.lnk.out_node, false, _innov.innovation_num2, 0);
            }
        }

        if (new_node == null) {
            //The innovation is totally novel
            //Create the new Genes
            //Create the new NNode
            //By convention, it will point to the first trait
            // get the current node id with postincrement
            int curnode_id = pop.getCur_node_id_and_increment();

            // pass this current nodeid to newnode and create the new node
            new_node = new NNode(NeatConstant.NEURON, curnode_id, NeatConstant.HIDDEN);
            new_node.nodetrait = traits.firstElement();

            // get the current gene inovation with post increment
            double gene_innov1 = pop.getCurr_innov_num_and_increment();

            // create gene with the current gene inovation
            newgene1 = new Gene(_gene.lnk.linktrait, 1.0, _gene.lnk.in_node, new_node, _gene.lnk.is_recurrent, gene_innov1, 0);

            // re-read the current innovation with increment
            double gene_innov2 = pop.getCurr_innov_num_and_increment();

            // create the second gene with this innovation incremented
            newgene2 = new Gene(_gene.lnk.linktrait, _gene.lnk.weight, new_node, _gene.lnk.out_node, false, gene_innov2, 0);

            pop.innovations.add(new Innovation(_gene.lnk.in_node.node_id, _gene.lnk.out_node.node_id, gene_innov1, gene_innov2, new_node.node_id, _gene.innovation_num));
        }

        //Now add the new NNode and new Genes to the Genome
        genes.add(newgene1);
        genes.add(newgene2);
        node_insert(nodes, new_node);
    }

    /**
     * Creation of a new random genome with :
     * new_id   = numerical identification of genome
     * i   = number of input nodes
     * o   = number of output nodes
     * n   = number of hidden nodes
     * nmax   = number max of node
     * this number must be >= (i + n + o)
     * r   = the network can have a nodes recurrent ?
     * linkprob = probability of a link from nodes ( must be in interval  ]0,1[);
     */
    public Genome(int new_id, int inputNodesCount, int outputNodesCount, int hiddenNodesCount, int nmax, boolean recurrentNodes, double linkprob) {
        double new_weight;

        NNode in_node = null;
        NNode out_node = null;

        notes = null;

        //
        //    i i i n n n n n n n n n n n n n n n n . . . . . . . . o o o o
        //    |                                   |                 ^     |
        //    |<----------- maxnode ------------->|                 |     |
        //    |                                                     |     |
        //    |<-----------------------total nodes -----------------|---->|
        //                                                          |
        //                                                          |
        //     first output ----------------------------------------+
        //
        //
        int totalnodes = inputNodesCount + outputNodesCount + nmax;

        traits = new Vector<>(Neat.p_num_trait_params, 0);
        nodes = new Vector<>(totalnodes, 0);
        genes = new Vector<>(totalnodes, 0);

        int matrixdim = totalnodes * totalnodes;

        boolean[] cm = new boolean[matrixdim]; //Dimension the connection matrix
        boolean[] cmp;

        int maxnode = inputNodesCount + hiddenNodesCount;
        int first_output = totalnodes - outputNodesCount + 1;

        //Assign the id
        genome_id = new_id;

        //Create a dummy trait (this is for future expansion of the system)
        Trait newtrait = new Trait(1, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        traits.add(newtrait);

        //Build the input nodes
        for (int count = 1; count <= inputNodesCount; count++) {
            NNode newnode = new NNode(NeatConstant.SENSOR, count, NeatConstant.BIAS);

            if (count < inputNodesCount) {
                newnode = new NNode(NeatConstant.SENSOR, count, NeatConstant.INPUT);
            }

            newnode.nodetrait = newtrait;
            //Add the node to the list of nodes
            nodes.add(newnode);
        }

        //Build the hidden nodes
        for (int count = inputNodesCount + 1; count <= inputNodesCount + hiddenNodesCount; count++) {
            NNode newnode = new NNode(NeatConstant.NEURON, count, NeatConstant.HIDDEN);
            newnode.nodetrait = newtrait;

            //Add the node to the list of nodes
            nodes.add(newnode);
        }

        //Build the output nodes
        for (int count = first_output; count <= totalnodes; count++) {
            NNode newnode = new NNode(NeatConstant.NEURON, count, NeatConstant.OUTPUT);
            newnode.nodetrait = newtrait;

            //Add the node to the list of nodes
            nodes.add(newnode);
        }

        boolean done = false;
        boolean rc1;
        boolean rc2;

        double forced_probability = 0.5;
        int abort = 0;

        while (!done) {
            abort++;
            if (abort >= 20) {
                linkprob = forced_probability;
                forced_probability += .01;
            }

            if (abort >= 700) {
                System.out.print("\n SEVERE ERROR in genome random creation costructor : genome has not created");
                System.exit(12);
            }

            //creation of connections matrix
            //Step through the connection matrix, randomly assigning bits
            cmp = cm;
            for (int count = 0; count < matrixdim; count++) {
                cmp[count] = NeatRoutine.randfloat() < linkprob;
            }

            //Connect the nodes
            int innov_number = 0; //counter for labelling the innov_num  of genes

            //Step through the connection matrix, creating connection genes
            cmp = cm;

            for (int col = 1; col <= totalnodes; col++) {
                for (int row = 1; row <= totalnodes; row++) {
                    if ((cmp[innov_number] && (col > inputNodesCount))
                            && ((col <= maxnode) || (col >= first_output))
                            && ((row <= maxnode) || (row >= first_output))) {
                        //If it isn't recurrent, create the connection no matter what

                        boolean create_gene = true;
                        boolean flag_recurrent = false;
                        if (col <= row) {
                            if (!recurrentNodes) {
                                create_gene = false;
                            }
                            flag_recurrent = true;
                        }

                        if (create_gene) {
                            int fnd = 0;

                            for (int i = 0; i < nodes.size() && fnd < 2; i++) {
                                NNode _node = nodes.get(i);
                                if (_node.node_id == row) {
                                    fnd++;
                                    in_node = _node;
                                }
                                if (_node.node_id == col) {
                                    fnd++;
                                    out_node = _node;
                                }
                            }

                            //Create the gene + link
                            new_weight = NeatRoutine.randposneg() * NeatRoutine.randfloat();
                            //Add the gene to the genome
                            genes.add(new Gene(newtrait, new_weight, in_node, out_node, flag_recurrent, innov_number, new_weight));
                        }
                    } // end condition for a correct link in genome
                    innov_number++;
                }
            }

            rc1 = verify();

            if (rc1) {
                Network net = this.genesis(genome_id);
                rc2 = net.is_minimal();

                if (rc2) {
                    int lx = net.max_depth();
                    int dx = net.is_stabilized(lx);

                    if (((dx == lx) && (!recurrentNodes)) || ((lx > 0) && (recurrentNodes) && (dx == 0))) {
                        done = true;
                    }
                }
                net.genotype = null;
                this.phenotype = null;
            }

            if (!done) {
                genes.clear();
            }
        }
    }

    public Genome(int id, IOseq xFile) {
        boolean done = false;

        genome_id = id;

        traits = new Vector<>(3, 0);
        nodes = new Vector<>(3, 0);
        genes = new Vector<>(3, 0);

        while (!done) {
            String xline = xFile.IOseqRead();
            StringTokenizer st = new StringTokenizer(xline);

            String curword = st.nextToken();
            if (curword.equalsIgnoreCase("genomeend")) {
                curword = st.nextToken();
                if (Integer.parseInt(curword) != genome_id) {
                    System.out.println(" *ERROR* id mismatch in genome");
                }
                done = true;
            } else if (curword.equals("/*")) {
                curword = st.nextToken();
                while (!curword.equals("*/")) {
                    curword = st.nextToken();
                }
            } else if (curword.equals("trait")) {
                Trait newtrait;
                newtrait = new Trait(xline);
                traits.addElement(newtrait);
            } else if (curword.equals("node")) {
                NNode newnode;
                newnode = new NNode(xline, traits);
                nodes.addElement(newnode);
            } else if (curword.equals("gene")) {
                Gene newgene;
                newgene = new Gene(xline, traits, nodes);
                genes.addElement(newgene);
            }
        }
    }

    public void print_to_file(IOseq xFile) {
        xFile.IOseqWrite("genomestart  " + genome_id);

        for (Trait _trait : traits) {
            _trait.print_to_file(xFile);
        }

        for (NNode _node : nodes) {
            _node.print_to_file(xFile);
        }

        for (Gene _gene : genes) {
            _gene.print_to_file(xFile);
        }

        xFile.IOseqWrite("genomeend " + genome_id);
    }

    public void View_mate_singlepoint(Genome g) {
        String mask4 = " 0000";
        DecimalFormat fmt4 = new DecimalFormat(mask4);

        int stopA;
        int stopB;
        int j;
        int j1;
        int j2;

        int crosspoint;

        Vector<Gene> genomeA;
        Vector<Gene> genomeB;

        if (genes.size() < g.genes.size()) {
            stopA = genes.size();
            stopB = g.genes.size();
            genomeA = genes;
            genomeB = g.genes;
        } else {
            stopA = g.genes.size();
            stopB = genes.size();
            genomeA = g.genes;
            genomeB = genes;
        }

        double[][] v3 = new double[g.genes.size() * 2][2];
        double[] vr = new double[g.genes.size() * 2];

        for (crosspoint = 0; crosspoint < stopA; crosspoint++) {
            int genecounter = 0; //Ready to count to crosspoint

            boolean done = false;
            double v1 = 0.0;
            double v2 = 0.0;
            j1 = 0;
            j2 = 0;
            j = 0;

            double cross_innov = 0;

            // compute what is the hight innovation
            double last_innovB = genomeB.elementAt(stopB - 1).innovation_num;

            while (!done) {
                boolean doneA = false;
                boolean doneB = false;

                if (j1 < stopA) {
                    v1 = genomeA.elementAt(j1).innovation_num;
                    doneA = true;
                }

                if (j2 < stopB) {
                    v2 = genomeB.elementAt(j2).innovation_num;
                    doneB = true;
                }

                if (doneA && doneB) {
                    if (v1 < v2) {
                        v3[j][0] = v1;
                        v3[j][1] = 0.0;
                        j1++;
                    } else if (v1 == v2) {
                        v3[j][0] = v1;
                        v3[j][1] = v1;
                        j1++;
                        j2++;
                    } else {
                        v3[j][0] = 0.0;
                        v3[j][1] = v2;
                        j2++;
                    }
                } else {
                    if (doneA) {
                        v3[j][0] = v1;
                        v3[j][1] = 0.0;
                        j1++;
                    } else if (doneB) {
                        v3[j][0] = 0.0;
                        v3[j][1] = v2;
                        j2++;
                    } else {
                        done = true;
                    }
                }

                if (!done) {
                    // innovA = innovB
                    if (v3[j][0] == v3[j][1]) {
                        if (genecounter < crosspoint) {
                            vr[j] = 1;
                            genecounter++;
                        } else if (genecounter == crosspoint) {
                            vr[j] = 3;
                            genecounter++;
                            cross_innov = v3[j][0];
                        } else if (genecounter > crosspoint) {
                            vr[j] = 2;
                            genecounter++;
                        }
                    } else if (v3[j][0] != 0 && v3[j][1] == 0) {
                        // innovA < innovB
                        if (genecounter < crosspoint) {
                            vr[j] = 1;                //  v3[j][0];
                            genecounter++;
                        } else if (genecounter == crosspoint) {
                            vr[j] = 1;                    // v3[j][1])
                            genecounter++;
                            cross_innov = v3[j][0];
                        } else if (genecounter > crosspoint) {
                            if (cross_innov > last_innovB) {
                                vr[j] = 1;
                                genecounter++;
                            }
                        }
                    } else if (v3[j][0] == 0 && v3[j][1] != 0) {
                        // innovA > innovB
                        if (genecounter < crosspoint) {
                            vr[j] = 0;                //  skip v3[j][0];
                        } else if (genecounter == crosspoint) {
                            vr[j] = 0;                    // skip
                        } else if (genecounter > crosspoint) {
                            if (cross_innov > last_innovB) {
                                vr[j] = 1;                    // v3[j][1];
                                genecounter++;
                            } else {
                                vr[j] = 2;
                                genecounter++;
                            }
                        }
                    }
                }

                j++;
            }

            int len_max = --j;

            // only for debug  : view innov's genomeA,B
            System.out.print("\n\n CROSSING SINGLE at index " + crosspoint);
            System.out.print("\n -- index -- ");
            int column = 0;
            for (j2 = 0; j2 < len_max; j2++) {
                if (v3[j2][0] > 0.0) {
                    System.out.print(fmt4.format(column++));
                } else {
                    System.out.print("     ");
                }
            }

            System.out.print("\n ----------- ");
            for (j2 = 0; j2 < len_max; j2++) {
                System.out.print("-----");
            }

            for (j1 = 0; j1 < 2; j1++) {
                System.out.print("\n Genome  [" + j1 + "] ");
                for (j2 = 0; j2 < len_max; j2++) {
                    System.out.print(fmt4.format((long) v3[j2][j1]));
                }
            }

            System.out.print("\n newgene [X] ");
            for (j2 = 0; j2 < len_max; j2++) {
                if (vr[j2] == 1) {
                    System.out.print("  AA ");
                } else if (vr[j2] == 2) {
                    System.out.print("  BB ");
                } else if (vr[j2] == 3) {
                    System.out.print("  XX ");
                } else if (vr[j2] == 4) {
                    System.out.print("  MM ");
                } else if (vr[j2] == 0) {
                    System.out.print("  -- ");
                }
            }

            System.out.print("\n");
        }
    }
}