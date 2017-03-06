
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#----------------------Nelder Mead Function----------------------------------------------------------------
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

NM_opt <- function( vtcs_init,            # initialiation of Nelder Mead vertices (to note that have to be in the range (-inf, +inf)
                    obj_fun,              # objective function
                    fxd_obj_param,        # objective function fixed parameter inputs
                    bdry_fun = NULL,      # boundary function to map vtcs for input into the selected algorithm
                    fxd_bdry_param = NULL,# parameters for input into the boundary functions
                    max_iter = 200,       # max iterations the algorithm will make in the search for local/global minima
                    max_0prgrss = 10,     # max iterations the algorithm will from the last improvement in objective value
                    vtcs_tol = 3,         # decimal places vertices tolerance before vtcs are considered equivalent
                    a_r = 1, a_e = 2, a_c = 0.5, a_s = 0.5 ) { #Nelder Mead parameters
                    #--------------------------------------
                    #Description of Nelder Mead parameters
                    #a_r -> reflect_cstnt
                    #a_e -> expand_cstnt
                    #a_c -> ctrt_cstnt
                    #a_s -> shrink_cstnt
                    
  #-------------------------------------------------------------------------------------
  # Nelder Mead specific function
  centroid <- function(all_vertices, worst_vertex) {
    output <- ( apply(all_vertices,2,sum) - worst_vertex ) / ( nrow(all_vertices) - 1 ) 
    output[1] <- round(output[1])
    output[2] <- round(output[2])
    return(output)
  }
  #-------------------------------------------------------------------------------------
  
  #checks if obj_fun is a function and passing 
  if (!is.function(obj_fun)) stop("input to obj_fun is not a function")
  
  #checks if unconstrained to constrained mapping is considered.
  if (!is.null(bdry_fun))  {
    boundary <- TRUE
    #checks if bdry_fun is a function
    if (!is.function(bdry_fun)) stop("input to bdry_fun is not a function")
    bdry_fun_args <- c(list(vtx = NULL),fxd_bdry_param)
  } else {
    boundary <- FALSE
  }
  
  #checks if obj_fun is a function and passing 
  if (!is.function(obj_fun)) stop("input to obj_fun is not a function")
  
  #checks if the vertices input are valid 
  if( !(is.matrix(vtcs_init) | is.data.frame(vtcs_init)) ) stop("vtcs_init format should be a matrix or data.frame")
  if( !length(dim(vtcs_init)) == 2 ) stop("dimensions of vtcs_init should be 4 vertices x 3 fields")
  if( !( (nrow(vtcs_init)- ncol(vtcs_init)) == 1) )  stop("[1] if optimization dimensions is N, number of vertices should be N + 1 
                                                           [2] please ensure that rows represent the vertices and columns the optization parameters")
  
  
  #passing in Nelder Mead vertices input and objtive function fixed parameters
  obj_fun_args <- c(list(vtx = NULL), fxd_obj_param)
  
  
  #passing in intialization parameters
  vtcs <- cbind( obj = rep(0,nrow(vtcs_init)), vtcs_init )
  
  #Miscellaneous
  N_vtcs <- nrow(vtcs)
  var_names <- colnames(vtcs)[!( colnames(vtcs) %in% c("obj") )]
  
  #initialization of vertices
  for ( i in 1:N_vtcs)  {
    if (boundary == TRUE) {
      bdry_fun_args[["vtx"]] <- vtcs[i, var_names]
      obj_fun_args[["vtx"]] <- do.call( bdry_fun, bdry_fun_args )
    }
    vtcs[i,"obj"] <- do.call( obj_fun, obj_fun_args )
  }
  vtcs <- vtcs[order(vtcs[,"obj"]),] # sort from best to worst
  
  
  #NM iteration initialization
  obj_tm1 <- obj_t <- 0
  N_shrink <- 0
  iter = 1
  iter_eq0 = 1
  
  
  #NM begins
	while ( (iter <= max_iter) & (iter_eq0 <= max_0prgrss) )	{
		path <- "reflection"
		X_w <- vtcs[N_vtcs,var_names] 
		C <- centroid(vtcs[,var_names], X_w ) #centroid of best hyperplane which is opposite the worst vertex 

		#---- Reflection step ----------------------------------------------------------------------------
		if (path == "reflection") {
		  X_r <- a_r*(C - X_w) 
		  if (boundary == TRUE) {
			bdry_fun_args[["vtx"]] <- X_r
			obj_fun_args[["vtx"]] <- do.call( bdry_fun, bdry_fun_args )
		  }
		  obj_r <- do.call( obj_fun, obj_fun_args )
		  #Relect if ( (vtcs[,"obj"] >= obj_r) & (coj_r > vtcs[,"obj"]) ) { 
		  if ( (vtcs[1,"obj"] <= obj_r) & (obj_r < vtcs[2,"obj"]) ) { 
			vtcs[N_vtcs,] <- c(cbj = obj_r, X_r); 
			# path <- "stopping" 
		  }	else if ( obj_r < vtcs[1,"obj"] ) { #go to expansion
			path <- "expansion" 
		  }	else {
			path <- "contraction"
		  }
		} # ---------------------------------------------------------------------------------------------


		#---- Expansion Step ----------------------------------------------------------------------------
		if (path == "expansion") {
		  X_e <- C + a_e*(X_r - C)
		  if (boundary == TRUE) {
			bdry_fun_args[["vtx"]] <- X_e
			obj_fun_args[["vtx"]] <- do.call( bdry_fun, bdry_fun_args )
		  }
		  obj_e <- do.call( obj_fun, obj_fun_args )
		  if (obj_e < obj_r)	{
			vtcs[N_vtcs,] <- c(obj = obj_e, X_e) 
		  } else { 
			vtcs[N_vtcs,] <- c(obj = obj_r, X_r)
		  }
		} # ---------------------------------------------------------------------------------------------


		#---- Contraction Step ---------------------------------------------------------------------------- 
		if (path == "contraction") { 
		  if ( (vtcs[2,"obj"] <= obj_r) & (obj_r < vtcs[N_vtcs,"obj"]) ) { # g(X_bad) <= g(X_r) < g(X_worst) 
			#outer contraction 
			X_o <- C + a_c*(X_r - C)
			if (boundary == TRUE) {
			  bdry_fun_args[["vtx"]] <- X_o
			  obj_fun_args[["vtx"]] <- do.call( bdry_fun, bdry_fun_args )
			}
			obj_o <- do.call( obj_fun, obj_fun_args )
			if (obj_o <= obj_r) { #g(X_o) <= g(X_r) 
			  vtcs[N_vtcs, ] <- c( obj = obj_o, X_o) 
			  #path <- "stopping" 
			} else { 
			  path <- "shrinking" 
			}
		  } else if (vtcs[N_vtcs,"obj"] <= obj_r) {# g(X_worst) <= g(X_r) 
			#inner contraction 
			X_i <- C + a_c*(X_w - C)
			if (boundary == TRUE) {
			  bdry_fun_args[["vtx"]] <- X_i
			  obj_fun_args[["vtx"]] <- do.call( bdry_fun, bdry_fun_args )
			}
			obj_i <- do.call( obj_fun, obj_fun_args )
			if (obj_i <= vtcs[N_vtcs,"obj"]) { #g(x_o) <= g(x_r) 
			  vtcs[N_vtcs, ] <- c( obj = obj_i, X_i) 
			  #path <- "stopping" 
			} else { 
			  path <- "shrinking" 
			} 
		  }
		} # ---------------------------------------------------------------------------------------------

		# Shrinking Step --------------------------------------------------------------------------------
		if (path == "shrinking") {
		  idx_s <- which(!vtcs[,"obj"] %in% vtcs[1,"obj"])
		  X_best <- t(matrix( rep(vtcs[1,var_names], length(idx_s)), nrow = N_vtcs - 1, ncol = length(idx_s) ))
		  vtcs[idx_s, var_names] <- X_best + a_s*(vtcs[idx_s,var_names] - X_best)
		  for (i in idx_s) {
			if (boundary == TRUE) {
			  bdry_fun_args[["vtx"]] <- vtcs[i, var_names]
			  obj_fun_args[["vtx"]] <- do.call( bdry_fun, bdry_fun_args )
			}
			vtcs[i,"obj"] <- do.call( obj_fun, obj_fun_args )
		  }
		}
		# -----------------------------------------------------------------------------------------------


		
		
		#---- Stopping Step -----------------------------------------------------------------------------
		vtcs <- vtcs[ order(vtcs[,"obj"]),] # reordering matrix according to obj
		obj_t <- vtcs[1,"obj"]
		
		# Mapping vertices back to orignal values for printing and returning output
		vtcs_print <- vtcs
		if (boundary == TRUE) {
      for ( i in 1:N_vtcs)  {
        bdry_fun_args[["vtx"]] <- vtcs[i, var_names]
        vtcs_print[i,var_names] <- do.call( bdry_fun, bdry_fun_args )
      }
		}
		
		flush.console()
		print('-----------------------------------------------------------')
		print(paste0('   Current iteration: ',iter))
		print('-----------------------------------------------------------')
		print(paste0('Objective delta: ', obj_t - obj_tm1) )
		print(paste0('Nelder Mead step taken: ', path))
		print('Current vertices coordinates and objective value')
		print(vtcs_print)
		print('-----------------------------------------------------------')
		print('')
		
		converged <- all( apply( vtcs_print[,-1], 2, function(x) var(round(x,2)) <= 0 ) )
		
		if (converged) {
		  iter_eq0 = max_0prgrss + 1
		} else if ( (obj_t - obj_tm1) != 0 )	{
		  iter_eq0 = 0
		}	else {
		  iter_eq0 = iter_eq0 + 1
		}        
		obj_tm1 <- obj_t
	  iter = iter + 1
	} 
  
  # Return results
  if (converged & (iter_eq0 > max_0prgrss)) {
    print('-----------------------------------------------------------')
    print('Algorithm converged to local/global minima.')
    print('-----------------------------------------------------------')
    
  } else if (!converged & (iter_eq0 > max_0prgrss)){
    print('-----------------------------------------------------------')
    print(paste0('Algorithm has ran for ', max_0prgrss, ' iterations without improvement'))
    print('-----------------------------------------------------------')
    
  } else if (iter < max_iter){
    print('-----------------------------------------------------------')
    print('Algorithm assumed to have converged.')
    print('-----------------------------------------------------------')
  } else {
    print('-----------------------------------------------------------')
    print(paste0('Algorithm exceeded specified max iterations of: ', max_iter))
    print('-----------------------------------------------------------')
  }
  
	return(vtcs_print)
}
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo    

















