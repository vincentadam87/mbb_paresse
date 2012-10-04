---- comments on the existing inversion scripts
---- October 2012  
invert_paresse_MATHYS: reward and effort learning with volatility 

invert_paresse_RAND: simple learning model without hidden states where three bias (hand, condition and session) are passed to the softmax decision rule. 

invert_paresse_RL: Q-learning reinforcement learning model and WSLS version (Q-learning model with the learning rate set to 1). Two versions of the discount function: linear and hyperbolic; 
		   two possibilities for the learning rate: shared between reward and effort learning or different. 

invert_paresse_STABLE_BAYES: the learning model that estimates posteriors based on the number of previous successful outcomes (big rewards and small efforts) for each hand. these posteriors are next passed to the softmax decision rule. Has 11 free parameters: beta, discount, bias + 8 parameters for the beta-distribution: (alpha,beta) x choice-dimension(reward-effort) x hand(left-right).

invert_paresse_mixed: models with different learnng rules for different choice dimensions (reward and effort). Learning rules could be: QL, WSLS, STABLE BAYES. Two discount function: linear and hyperbolic. For the model where both learning rules are QL (similar to invert_paresse_RL): test the possibility for shared or different learning rates. No bias included.    


