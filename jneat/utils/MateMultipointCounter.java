package jneat.utils;

import jNeatCommon.NeatRoutine;
import jneat.Gene;

import java.util.Vector;

public class MateMultipointCounter {
    boolean disable = false;
    boolean skip = false;
    boolean firstGenomeBetter = false;

    Gene currentGene = null;
    Gene chosenGene = null;

    int firstGenePosition = 0;
    int secondGenePosition = 0;

    Vector<Gene> firstGenes;
    Vector<Gene> secondGenes;

    public MateMultipointCounter(boolean firstGenomeBetter, Vector<Gene> firstGenes, Vector<Gene> secondGenes) {
        this.firstGenomeBetter = firstGenomeBetter;
        this.firstGenes = firstGenes;
        this.secondGenes = secondGenes;
    }

    public void reset() {
        skip = false;
    }

    public int getFirstGenePosition() {
        return firstGenePosition;
    }

    public int getSecondGenePosition() {
        return secondGenePosition;
    }

    public Gene getChosenGene() {
        return chosenGene;
    }

    public boolean shouldSkip() {
        return skip;
    }

    public boolean shouldDisable() {
        return disable;
    }

    public void resetDisable() {
        disable = false;
    }

    public void chooseGeneFromFirst() {
        chosenGene = firstGenes.elementAt(firstGenePosition);
        firstGenePosition++;
        if (!firstGenomeBetter) {
            skip = true; //Skip excess from the worse genome
        }
    }

    public void chooseGeneFromSecond() {
        chosenGene = secondGenes.elementAt(secondGenePosition);
        secondGenePosition++;
        if (firstGenomeBetter) {
            skip = true; //Skip excess from the worse genome
        }
    }

    public void chooseBetterGene() {
        Gene firstGene = firstGenes.elementAt(firstGenePosition);
        Gene secondGene = secondGenes.elementAt(secondGenePosition);

        if (firstGene.getInnovation_num() == secondGene.getInnovation_num()) {
            chosenGene = NeatRoutine.randfloat() < 0.5 ? firstGene : secondGene;

            // If one is disabled, the corresponding gene in the offspring will likely be disabled
            disable = (!firstGene.getEnable() || !secondGene.getEnable()) && NeatRoutine.randfloat() < 0.75;
            firstGenePosition++;
            secondGenePosition++;
        } else if (firstGene.getInnovation_num() < secondGene.getInnovation_num()) {
            chosenGene = firstGene;
            firstGenePosition++;
            if (!firstGenomeBetter) {
                skip = true;
            }
        } else {
            chosenGene = secondGene;
            secondGenePosition++;
            if (firstGenomeBetter) {
                skip = true;
            }
        }
    }

    public void checkDuplicate(Vector<Gene> newGenes) {
        skip |= newGenes.stream().anyMatch(
                gene -> gene.getLnk().getIn_node().getNode_id() == chosenGene.getLnk().getIn_node().getNode_id()
                && gene.getLnk().getOut_node().getNode_id() == chosenGene.getLnk().getOut_node().getNode_id()
                && gene.getLnk().getIs_recurrent() == chosenGene.getLnk().getIs_recurrent()

                || gene.getLnk().getIn_node().getNode_id() == chosenGene.getLnk().getOut_node().getNode_id()
                && gene.getLnk().getOut_node().getNode_id() == chosenGene.getLnk().getIn_node().getNode_id()
                && !gene.getLnk().getIs_recurrent()
                && !chosenGene.getLnk().getIs_recurrent()
        );
    }
}
