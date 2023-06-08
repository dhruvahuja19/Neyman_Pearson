from decimal import Decimal
import math
import argparse
import matplotlib.pyplot as plt

def disc_NP(Dist1, Dist2, alpha):
    """
    This function takes two distributions and a significance level and returns the decision rule for the Neyman-Pearson test.
    Args:
    Dist1: The first distribution as a list of tuples (value, probability)
    Dist2: The second distribution as a list of tuples
    alpha: The significance level
    :var map1: A dictionary with the values of Dist1 as keys and the probabilities of Dist1 as values(x, P(X = x))
    :var map2: A dictionary with the values of Dist2 as keys and the probabilities of Dist2 as values(x, P(X = x))
    :var likelihood_ratio: A list of tuples (y, L(y)) where y a possible value and L(y) is the likelihood ratio
    :var map3: A dictionary with potential thresholds as keys and probabilities of those thresholds as values
    under the null hypothesis (t, P(L(y) = t | H0))
    :var map4: A dictionary with potential thresholds as keys and probabilities of those thresholds as values under
    the alt hyptohsis (t, P(L(y) = t | H1))
    :var threshold: self explan
    :var Ly_to_y: A dictionary with likelihood ratios as keys and the corresponding values of y as values(Used for printing results)
    """
    def create_points(map3, map4, gamma):
        """This function will return a list of triples (x, y, z) where x is the threshold,
         y is the probability of false alarm, and z is the probability of detection
         """
        # TO calculate the PFA of some threshold t, then we need to calculate P(L(y) > t | H0) which
        # is the sum of map3 values from t to the end of the list
        # For PCD, we care about P(L(y) > t | H1) which requires the declaration of a new map4. Map3 is based on Dist1
        # and map4 is based on Dist2. We can use the same keys for both maps, but the values will be different.
        points = []
        assert len(map3) == len(map4)
        keys = list(map3.keys())
        keys.sort()
        for i in range(0, len(keys)):
            t = keys[i]
            pfa = Decimal(sum(list(map3.values())[i + 1:])) + Decimal(gamma) * Decimal(map3[t])
            pcd = Decimal(sum(list(map4.values())[i + 1:])) + Decimal(gamma) * Decimal(map4[t])
            points.append((pfa, pcd, t))
        return points
    def print_results(threshold, gamma, nullvals, altvals, gammavals):
        print("Accept null hypothesis if L(y) <", threshold)
        print("Reject null hypothesis if L(y) >", threshold)
        print("If L(y) = ", threshold, "then reject the null hypothesis with probability", gamma)
        print("More formally, the decision rule is: If y is in",
              nullvals, "then accept the null hypothesis")
        print("If y is in", altvals,
              "then reject the null hypothesis")
        print("If y is in", gammavals,
              "then reject the null hypothesis with probability", gamma)
        return "The decision rule is: If y is in: " + str(nullvals) + " then accept the null hypothesis" + "\n" + \
                    "If y is in: " + str(altvals) + " then reject the null hypothesis" + "\n" + \
                    "If y is in: " + str(gammavals) + " then reject the null hypothesis with probability: " + str(gamma)
    dist1_range = [val for val, prob in Dist1]
    dist2_range = [val for val, prob in Dist2]
    map1, map2, likelihood_ratio = {}, {}, []
    for val, prob in Dist1:
        map1[val] = prob
    for val, prob in Dist2:
        map2[val] = prob
    [Dist2.append((val, 0)) for val in dist1_range if val not in dist2_range]
    [Dist1.append((val, 0)) for val in dist2_range if val not in dist1_range]
    assert len(Dist1) == len(Dist2)
    for val, prob in Dist1:
        if map2.get(val) == 0:
            likelihood_ratio.append((val, math.inf))
        else:
            likelihood_ratio.append((val, round(Decimal(map2.get(val)) / Decimal(prob), 6)))
    likelihood_ratio.sort(key=lambda x: x[1])
    map3, map4 = {}, {}
    Ly_to_y = {}
    for y, Ly in likelihood_ratio:
        if Ly_to_y.get(Ly) == None:
            Ly_to_y[Ly] = [y]
        else:
            Ly_to_y[Ly].append(y)
    for i in range(len(likelihood_ratio)):
        if map3.get(likelihood_ratio[i][1]) == None:
            map3[likelihood_ratio[i][1]] = map1.get(likelihood_ratio[i][0])
        else:
            map3[likelihood_ratio[i][1]] = round(Decimal(map3.get(likelihood_ratio[i][1]))
                                                 + Decimal(map1.get(likelihood_ratio[i][0])), 6)
    for i in range(len(likelihood_ratio)):
        if map4.get(likelihood_ratio[i][1]) == None:
            map4[likelihood_ratio[i][1]] = map2.get(likelihood_ratio[i][0])
        else:
            map4[likelihood_ratio[i][1]] = round(Decimal(map4.get(likelihood_ratio[i][1]))
                                                 + Decimal(map2.get(likelihood_ratio[i][0])), 6)
    i = 1
    keys = list(map3.keys())
    keys.sort()
    curr_threshold = likelihood_ratio[0][1]
    curr_pfa = round(sum([round(Decimal(map1.get(likelihood_ratio[i][0])), 6) for i in range(1, len(likelihood_ratio))
                    if likelihood_ratio[i][1] > likelihood_ratio[0][1]]), 6)
    while curr_pfa > alpha and i < len(keys):
        curr_threshold = keys[i]
        curr_pfa -= round(Decimal(map3.get(keys[i])), 6)
        i += 1
    nullvals = [[y for y in Ly_to_y[Ly]] for Ly in list(Ly_to_y.keys()) if Ly < curr_threshold]
    altvals = [[y for y in Ly_to_y[Ly]] for Ly in list(Ly_to_y.keys()) if Ly > curr_threshold]
    gammavals = [[y for y in Ly_to_y[Ly]] for Ly in list(Ly_to_y.keys()) if Ly == curr_threshold]
    if curr_pfa == alpha:
        answer = print_results(curr_threshold, 0, nullvals, altvals, gammavals)
    else:
        gamma = round(round(Decimal(alpha) - Decimal(curr_pfa), 6) / round(Decimal(map3.get(curr_threshold)), 6), 6)
        answer = print_results(curr_threshold, gamma, nullvals, altvals, gammavals)
    points = create_points(map3, map4, gamma)
    plot_roc_curve_with_answer(points, answer)
    return [(0, nullvals), (1, altvals), (gamma, gammavals)]
def Continous_Neyman_Pearson(Dist1, Dist2, alpha):
    """
    This function takes two distributions and a significance level and returns the decision rule for the Neyman-Pearson test.
    Args:
    Dist1: The first distribution(in the form of a function, the pdf)
    Dist2: The second distribution(in the form of a function, the pdf)
    alpha: The significance level
    """
    # Instead of calculating L(y) for all y which is impossible due to continuity
    # I can instead divide the pdfs to create L(y)


def plot_roc_curve_with_answer(points, answer):
    # Ensure the input list is not empty
    if not points:
        raise ValueError("Input list cannot be empty.")

    # Separate the PFA, PCD, and labels from the input points
    pfas, pcds, labels = zip(*points)

    # Add (0, 0) and (1, 1) points to represent the lower-left and upper-right corners
    pfas = (0,) + pfas + (1,)
    pcds = (0,) + pcds + (1,)
    labels = ["(0, 0)"] + list(labels) + ["(1, 1)"]

    # Sort the points based on PFA in ascending order
    sorted_indices = sorted(range(len(pfas)), key=lambda i: pfas[i])
    pfas = [pfas[i] for i in sorted_indices]
    pcds = [pcds[i] for i in sorted_indices]
    labels = [labels[i] for i in sorted_indices]

    # Calculate the area under the ROC curve (AUC)
    auc = sum((pfas[i + 1] - pfas[i]) * pcds[i] for i in range(len(pfas) - 1))

    # Plot the ROC curve and add labels
    plt.plot(pfas, pcds, marker='o', linestyle='-')
    plt.xlabel('PFA (Probability of False Alarm)')
    plt.ylabel('PCD (Probability of Correct Detection)')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.grid(True)
    plt.text(0.5, 0.2, f'AUC = {auc:.2f}', ha='center', fontsize=12)

    for i, label in enumerate(labels):
        plt.annotate(label, (pfas[i], pcds[i]), textcoords="offset points", xytext=(5,5), ha='center')

    # Create a textbox to display the string "answer"
    plt.subplots_adjust(bottom=0.2)
    textbox = plt.axes([0.15, 0.05, 0.7, 0.1])
    textbox.text(0.5, 0.5, answer, fontsize=6, ha='center', va='center')

    plt.show()

def main ():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--Null',
        help='The distribution of the null hypothesis as a list of tuples [(value, P[Value])...].',
        type= parse_tuple_list,
        required=False)
    parser.add_argument(
        '--Alt',
        help='The distribution of the hypothesis as a list of tuples [(value, P[Value])...].',
        type=parse_tuple_list,
        required=False)
    parser.add_argument(
        '--alpha', help='Alpha', type=float, required=False)
    args = parser.parse_args()
    disc_NP(args.Null_hypothesis_distribution, args.Alternative_hypothesis_distribution, args.alpha)

def parse_tuple_list(arg):
    try:
        tuples = eval(arg)  # Safely evaluate the argument as a Python expression
        if isinstance(tuples, list) and all(isinstance(t, tuple) and len(t) == 2 for t in tuples):
            return tuples
        else:
            raise argparse.ArgumentTypeError("Invalid tuple list format.")
    except:
        raise argparse.ArgumentTypeError("Invalid tuple list format.")

if __name__ == '__main__':
  main()



