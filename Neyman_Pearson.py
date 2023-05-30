from decimal import Decimal
def disc_NP(Dist1, Dist2, alpha):
    """
    This function takes two distributions and a significance level and returns the decision rule for the Neyman-Pearson test.
    Args:
    Dist1: The first distribution as a list of tuples (value, probability)
    Dist2: The second distribution as a list of tuples
    alpha: The significance level
    """
    def print_results(threshold, gamma, Ly_to_y, nullvals, altvals, gammavals):
        print("Accept null hypothesis if L(y) <", threshold)
        print("Reject null hypothesis if L(y) >", threshold)
        print("If L(y) = ", threshold, "then reject the null hypothesis with probability", gamma)
        print("More formally, the decision rule is: If y is in",
              nullvals, "then accept the null hypothesis")
        print("If y is in", altvals,
              "then reject the null hypothesis")
        print("If y is in", gammavals,
              "then reject the null hypothesis with probability", gamma)

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
        likelihood_ratio.append((val, map2.get(val) / prob))
    likelihood_ratio.sort(key=lambda x: x[1])
    map3 = {}
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
            map3[likelihood_ratio[i][1]] = map3.get(likelihood_ratio[i][1]) + map1.get(likelihood_ratio[i][0])
    i = 1
    keys = list(map3.keys())
    keys.sort()
    curr_threshold = likelihood_ratio[0][1]
    curr_pfa = sum([map1.get(likelihood_ratio[i][0]) for i in range(1, len(likelihood_ratio))
                    if likelihood_ratio[i][1] > likelihood_ratio[0][1]])
    while curr_pfa > alpha and i < len(keys):
        curr_threshold = keys[i]
        curr_pfa -= map3.get(keys[i])
        i += 1
    nullvals = [[y for y in Ly_to_y[Ly]] for Ly in list(Ly_to_y.keys()) if Ly < curr_threshold]
    altvals = [[y for y in Ly_to_y[Ly]] for Ly in list(Ly_to_y.keys()) if Ly > curr_threshold]
    gammavals = [[y for y in Ly_to_y[Ly]] for Ly in list(Ly_to_y.keys()) if Ly == curr_threshold]
    if curr_pfa == alpha:
        print_results(curr_threshold, 0, Ly_to_y, nullvals, altvals, gammavals)
    else:
        gamma = (alpha - curr_pfa) / map3.get(curr_threshold)
        print_results(curr_threshold, gamma, Ly_to_y, nullvals, altvals, gammavals)
    return [(0, nullvals), (1, altvals), (gamma, gammavals)]


