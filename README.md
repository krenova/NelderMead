# NelderMead

A Generic Nelder Mead optimization function, _**NM_opt()**_, for the parameter optimization of any predictive algorithms that is or can be integrated with R. This means that unlike for example, R's **caret** grid search function, the algorithm of choice does not need to be integrated or recognized by **caret** but as long as the algorithm can run in R,   _**NM_opt()**_, can be used to find the algorithm's local optimal parameters. The extent of this means that even python algorithms can be optimized using  _**NM_opt()**_. However, in exchange for this flexibility,  _**NM_opt()**_ is not easy to use and will require sometime to understand how the function is built and how it can be fully exploited.


1. [NelderMead](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) is a heuristic unconstrained optimization algorithm that is loosely speaking, based around the construction of a simplex (the simplest plane that be constructed in higher dimensions), that is used to get a sense of the optimization surface. Based on some logical rules, this simplex will then 'feel' and work its way towards a local minima/maxima.


2. To note that because NelderMead is an unconstrained optimization algorithm, mapping and inverse mapping of parameters that are constrained will be required. For example, parameters restricted to lie in the range (0, +inf) will need a log transformation to map the range to the values (-inf, inf).
 
3. For a brief overview of the function and the algorithm's origins, please refer to the following document
	* [document.pdf](https://github.com/krenova/NelderMead/blob/master/document.pdf)
 	* to note that the document is slightly outdated and may not reflect all the current function's iterations
 	* for the updated function arguments, please refer to the code itself: [NelderMead.R](https://github.com/krenova/NelderMead/blob/master/lib/NelderMead.R)

 
4. Before attempting to use the function, going through either the following examples is highly recommended:
 	* to quickly jump into the use of the library, please refer to: [example.ipynb](https://github.com/krenova/NelderMead/blob/master/example.ipynb)
 	* to run the programme in R instead of Jupyter Notebook, please refer to: [example.R](https://github.com/krenova/NelderMead/blob/master/example.R)

5. _**NM_opt()**_ is built on base R and therefore only has dependencies on changes in base R code. However, to note that in the examples provided where _**NM_opt()**_ was applied on **xgboost** and **randomForest** methods, changes in their functions may cause the examples to be unable to run. In this case, to run the codes, either modify the respective wrappers around their functions accordingly or revert the packages to versions in 2017-03-01.