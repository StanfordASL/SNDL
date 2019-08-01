function [vars,prob] = Newton_descent(vars_prev,N_steps,term_tol,line_alpha,line_beta,t_lb,obj,JH_obj,eps_H,prob)


%%
%initialize
it = 0;
converged = 0;
stalled = 0;

%% Newton descent

while (it<N_steps) && (~converged) && (~stalled)
    
    %% Get Descent direction
    
    [obj_val,prob] = obj(vars_prev,prob);
    [grad,Hobj] = JH_obj(vars_prev,prob);
    
    if (prob.check_W && max(prob.MAX_W) < 0.0)
        disp('Feasible!');
        vars = vars_prev;
        return;
    end
        
    step = -(Hobj + eps_H)\grad;

 %% Backtracking line search
 
    dec = - grad'*step;
    
    if (dec/2 > term_tol)

        step_accept = 0;
        t = 1;
        
        while (~step_accept) && (t>t_lb)
            
            vars_t = vars_prev + t*step;
            
            [obj_t,prob] = obj(vars_t,prob);
            
            %check cost decrease
            if (obj_t <= obj_val + line_alpha*t*(grad'*step))
                
                %cost decrease good
                step_accept = 1;
                
               if (prob.check_W) 
                   %only for metric infeas problem
                   if (max(prob.MAX_W) < 0.0)
                        %terminate
                        disp('Feasible!');
                        fprintf('it: %d, obj: %.4f, max F viol: %.4f, max W viol: %.4f, |grad|: %.4f, diff: %.6f, step: %.6f \n', ...
                            it+1, obj_t,max(prob.MAX_F),max(prob.MAX_W),norm(grad), norm(vars_t-vars_prev,'inf'),t);
                        
                        %save and return
                        vars = vars_t;
                        return;
                   end
               end
            else
                t = line_beta*t;
            end
        end
    else
        %descent no longer possible -> converged
        converged = 1;
    end
        
    %% Take step
    
    if (converged) 
        vars = vars_prev;
        diff = 0;
        t = 0;
    else    
        if (step_accept)
            
            vars = vars_t; %new set of values good
            obj_val = obj_t;
            diff = norm(vars-vars_prev,'inf');
            vars_prev = vars;
            
            if diff < term_tol
                converged = 1;
            end
        else
            %backtrack failed -> stalled
            vars = vars_prev;
            [~,prob] = obj(vars,prob);
            diff = 0;
            stalled = 1;
        end
    end
    
    %% Print

    it = it+1;
    fprintf('it: %d, obj: %.4f, max F viol: %.4f, max W viol: %.4f, |grad|: %.4f, diff: %.6f, step: %.6f \n', ...
            it, obj_val,max(prob.MAX_F),max(prob.MAX_W),norm(grad), diff,t);
    
end

%% Final print

if (converged)
    fprintf('Converged!\n');
elseif (it==N_steps)
    fprintf('Max iterations reached!\n');
else
    fprintf('Stalled!\n');
end

end