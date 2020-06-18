--- ImpactModeling ---

Andy, Joah, Natalie and Phil advised by Will 
Dr. Posa 
DAIR Lab

--- Recreating Models ---

Wang Mason Model:
Our first goal was to understand and create a function that applies the Wang Mason model for impacts
It takes in a given pre-impact state and returns the post impact state.

We then took our code for the Wang Mason model and applied it to data from Nima Fazeli, trying to find the optimal
mu and epsilon values to minimize the model's error.
We looked at both dropping a square (s = 0.06m) as well as a an ellipse (Maj Axis = 0.07m, Minor Axis = 0.05m)

After this, we compared our function with Nima's since our apporaches were different (comapre Wang Folder).

 
Anitescu Potra Poisson (AP Poisson) Model:
Next, we created a similar function with similar goals as the Wang function, but instead uses the AP Poisson model. 

LCP Solver Used:
Andreas Almqvist (2020). A pivoting algorithm solving linear complementarity problems 
(https://www.mathworks.com/matlabcentral/fileexchange/41485-a-pivoting-algorithm-solving-linear-complementarity-problems), 
MATLAB Central File Exchange. Retrieved June 18, 2020.

Again, we compared our function with Nima's and optimized to find the best mu and epsilon.

Ideal Rigid Body Bound (IRB Bound):

Unlike the prior two models, this model is not a predictive model. This model optimizes over the impulse space that is 
energy alloweable, and picks the impulses which lead to the lowest error comapred to the data. 
Using fmincon and non linear constraints, we wrote code that executed the IRB model.

--- Improving Models ---

Our next goal was to look at what assumptions could be improved and make impact models more accurate.

Point Contact:
First, we considered the case where the impact does not occour at a single point, however maybe over a patch. 
This then allows for a moment to be created by the shifting the location of the normal force. 




 


