--- /tmp/global_best.py	2023-03-10 00:31:21.212185550 +0530
+++ /home/sandeep/Software/anaconda3/lib/python3.7/site-packages/pyswarms/single/global_best.py	2023-03-10 00:30:30.692159350 +0530
@@ -213,7 +213,7 @@
             self.swarm.best_pos, self.swarm.best_cost = self.top.compute_gbest(self.swarm)
             # fmt: on
             if verbose:
-                self.rep.hook(best_cost=self.swarm.best_cost)
+                self.rep.hook(best_cost=str(self.swarm.best_cost))
             # Save to history
             hist = self.ToHistory(
                 best_cost=self.swarm.best_cost,
@@ -227,7 +227,7 @@
             relative_measure = self.ftol * (1 + np.abs(best_cost_yet_found))
             delta = (
                 np.abs(self.swarm.best_cost - best_cost_yet_found)
-                < relative_measure
+                < self.ftol
             )
             if i < self.ftol_iter:
                 ftol_history.append(delta)
