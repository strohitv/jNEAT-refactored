package jneat.utils;

public class CompabilityCounter {
    private int firstGenePosition = 0;
    private int secondGenePosition = 0;

    private int matchingCount = 0;
    private int disjointCount = 0;
    private int excessCount = 0;

    private double totalWeigthDifference = 0.0;

    public void addOwnExcess() {
        excessCount += 1;
        firstGenePosition++;
    }

    public void addOtherExcess() {
        excessCount += 1;
        secondGenePosition++;
    }

    public void addMatching(double firstMutationNumber, double secondMutationNumber) {
        matchingCount += 1.0;
        totalWeigthDifference += Math.abs(firstMutationNumber - secondMutationNumber);
        firstGenePosition++;
        secondGenePosition++;
    }

    public void addOwnDisjoint() {
        disjointCount += 1;
        firstGenePosition++;
    }

    public void addOtherDisjoint() {
        disjointCount += 1;
        secondGenePosition++;
    }

    public int getFirstGenePosition() {
        return firstGenePosition;
    }

    public int getSecondGenePosition() {
        return secondGenePosition;
    }

    public double getDisjointCount() {
        return disjointCount;
    }

    public double getExcessCount() {
        return excessCount;
    }

    public double getMatchingCount() {
        return matchingCount;
    }

    public double getTotalWeigthDifference() {
        return totalWeigthDifference;
    }
}
