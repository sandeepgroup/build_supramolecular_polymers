--- /home/ap/Arunima/anaconda3/lib/python3.9/site-packages/pyswarms/single/global_best.py	2023-03-18 00:58:43.284464768 +0530
+++ Global_best_updated.py	2023-03-18 00:59:01.388921267 +0530
@@ -204,6 +204,8 @@
 
         self.swarm.pbest_cost = np.full(self.swarm_size[0], np.inf)
         ftol_history = deque(maxlen=self.ftol_iter)
+        cost_history=[]
+        pos_history=[]
         #global cost_history
         for i in self.rep.pbar(iters, self.name) if verbose else range(iters):
             # Compute cost for current position and personal best
@@ -224,6 +226,8 @@
                 position=self.swarm.position,
                 velocity=self.swarm.velocity,
             )
+            cost_history.append(self.swarm.best_cost)
+            pos_history.append(self.swarm.best_pos.tolist())
             self._populate_history(hist)
             # Verify stop criteria based on the relative acceptable cost ftol
             relative_measure = self.ftol * (1 + np.abs(best_cost_yet_found))
