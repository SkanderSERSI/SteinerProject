package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

public class DefaultTeam {

	public ArrayList<Object> calculShortestPaths(ArrayList<Point> points, int edgeThreshold) {

		/** Exercice 1 **/
		ArrayList<Object> res = new ArrayList<Object>();
		int[][] paths = new int[points.size()][points.size()];
		double[][] dist = new double[points.size()][points.size()];

		for (int i = 0; i < dist.length; i++) {
			for (int j = 0; j < dist.length; j++) {
				if (points.get(i).distance(points.get(j)) < edgeThreshold) {
					dist[i][j] = points.get(i).distance(points.get(j));
					paths[i][j] = j;
				} else {
					dist[i][j] = Double.POSITIVE_INFINITY;
					paths[i][j] = -1;
				}
			}
		}

		for (int k = 0; k < dist.length; k++) {
			for (int i = 0; i < dist.length; i++) {
				for (int j = 0; j < dist.length; j++) {
					if (dist[i][j] > dist[i][k] + dist[k][j]) {
						dist[i][j] = dist[i][k] + dist[k][j];
						paths[i][j] = paths[i][k];
					}
				}
			}
		}
		res.add(paths);
		res.add(dist);
		return res;
	}

	public ArrayList<Edge> kruskal(ArrayList<Point> hitPoints, double[][] dist,
			HashMap<Point, Integer> mapPointIndice) {

		ArrayList<Edge> edges = new ArrayList<Edge>();
		for (Point p : hitPoints) {
			for (Point q : hitPoints) {
				if (p.equals(q) || contains(edges, p, q))
					continue;
				edges.add(new Edge(p, q));
			}
		}
		edges = sort(edges, dist, mapPointIndice);
		ArrayList<Edge> kruskal = new ArrayList<Edge>();
		Edge current;
		NameTag forest = new NameTag(hitPoints);
		while (edges.size() != 0) {
			current = edges.remove(0);
			if (forest.tag(current.p) != forest.tag(current.q)) {
				kruskal.add(current);
				forest.reTag(forest.tag(current.p), forest.tag(current.q));
			}
		}

		return kruskal;
	}

	private boolean contains(ArrayList<Edge> edges, Point p, Point q) {
		for (Edge e : edges) {
			if (e.p.equals(p) && e.q.equals(q) || e.p.equals(q) && e.q.equals(p))
				return true;
		}
		return false;
	}

	@SuppressWarnings("unchecked")
	private Tree2D edgesToTree(ArrayList<Edge> edges, Point root) {
		ArrayList<Edge> remainder = new ArrayList<Edge>();
		ArrayList<Point> subTreeRoots = new ArrayList<Point>();
		Edge current;
		while (edges.size() != 0) {
			current = edges.remove(0);
			if (current.p.equals(root)) {
				subTreeRoots.add(current.q);
			} else {
				if (current.q.equals(root)) {
					subTreeRoots.add(current.p);
				} else {
					remainder.add(current);
				}
			}
		}

		ArrayList<Tree2D> subTrees = new ArrayList<Tree2D>();
		for (Point subTreeRoot : subTreeRoots)
			subTrees.add(edgesToTree((ArrayList<Edge>) remainder.clone(), subTreeRoot));

		return new Tree2D(root, subTrees);
	}

	private ArrayList<Edge> sort(ArrayList<Edge> edges, double[][] dist, HashMap<Point, Integer> map) {

		if (edges.size() == 1) {

			return edges;
		}

		ArrayList<Edge> left = new ArrayList<Edge>();
		ArrayList<Edge> right = new ArrayList<Edge>();
		int n = edges.size();
		for (int i = 0; i < n / 2; i++) {
			left.add(edges.remove(0));
		}
		while (edges.size() != 0) {
			right.add(edges.remove(0));
		}
		left = sort(left, dist, map);
		right = sort(right, dist, map);

		ArrayList<Edge> result = new ArrayList<Edge>();
		while (left.size() != 0 || right.size() != 0) {
			if (left.size() == 0) {
				result.add(right.remove(0));
				continue;
			}
			if (right.size() == 0) {
				result.add(left.remove(0));
				continue;
			}
			if (dist[map.get(left.get(0).p)][map.get(left.get(0).q)] < dist[map.get(right.get(0).p)][map
					.get(right.get(0).q)])
				result.add(left.remove(0));
			else
				result.add(right.remove(0));
		}
		return result;
	}

	public Point getKeyFromValue(HashMap<Point, Integer> map, int value) {
		for (Entry<Point, Integer> entry : map.entrySet()) {
			if (entry.getValue() == value)
				return entry.getKey();
		}
		// System.out.println("hahahahahhahahahahahhahahaha");
		return new Point(-1, -1);
	}

	public Tree2D calculSteiner(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
		ArrayList<Object> path_dist = calculShortestPaths(points, edgeThreshold);
		int[][] paths = (int[][]) path_dist.get(0);
		double[][] dist = (double[][]) path_dist.get(1);

		HashMap<Point, Integer> mapPointIndice = new HashMap<>();

		for (int i = 0; i < points.size(); i++) {
			mapPointIndice.put(points.get(i), i);
		}

		ArrayList<Edge> kruskhaled = kruskal(hitPoints, dist, mapPointIndice);

		ArrayList<Edge> best_edges = new ArrayList<>();
		for (Edge edge : kruskhaled) {
			Point p = edge.p;
			Point q = edge.q;

			ArrayList<Edge> edges_intermediaires = new ArrayList<>();

			int i = mapPointIndice.get(p);
			int j = mapPointIndice.get(q);

			if (paths[i][j] != j) {

				int k = paths[i][j];
				edges_intermediaires
						.add(new Edge(getKeyFromValue(mapPointIndice, i), getKeyFromValue(mapPointIndice, k)));
				while (k != j) {
					int k_moins = k;
					k = paths[k][j];
					edges_intermediaires.add(
							new Edge(getKeyFromValue(mapPointIndice, k_moins), getKeyFromValue(mapPointIndice, k)));

				}
				edges_intermediaires
						.add(new Edge(getKeyFromValue(mapPointIndice, k), getKeyFromValue(mapPointIndice, j)));
				best_edges.addAll(edges_intermediaires);

			} else {
				best_edges.add(edge);
			}

		}
		System.out.println(best_edges);
		return edgesToTree(best_edges, best_edges.get(0).p);

	}
	
	public int compter(ArrayList<Point> points,Point p) {
		int cpt = 0;
		for (Point point : points) {
			if (point.x == p.x && point.y == p.y) {
				//System.out.println("point.x = "+point.x+" point.y = "+point.y);
				//System.out.println("p.x = "+p.x+" p.y = "+p.y);
				cpt++;
			}
		}
		return cpt;
		
	}
	
	// Probleme il y a des points des edges qui existent plus que prevu 
	// 1 edge d'extremite a au moins un de ses points qui existe qu'une seulee fois dans 
	// l'array list eclatee
	public ArrayList<Edge> filtrage_edges_extrems(ArrayList<Edge> edges){
		ArrayList<Edge> edgesAlgo = (ArrayList<Edge>) edges.clone();
		ArrayList<Point> points = eclatement_edges(edges);
		/*for (Point point : points) {
			System.out.println("("+point.x+","+point.y+")");
		}*/
		for (Edge edge1 : edges) {

			int cp1 = compter(points, edge1.p);
			int cp2 = compter(points, edge1.q);
			
			if (cp1 <= 1 || cp2 <= 1) {
				System.out.println("cp1 = " + cp1 + " cp2= " + cp2);
				edgesAlgo.remove(edge1);
				break;
			}
		}
		return edgesAlgo;	
	}
	public ArrayList<Point> filtrage_sommets_extrems(ArrayList<Edge> edges){
		ArrayList<Edge> edgesAlgo = (ArrayList<Edge>) edges.clone();
		ArrayList<Point> points = eclatement_edges(edges);
		ArrayList<Point> res = new ArrayList<>();
		/*for (Point point : points) {
			System.out.println("("+point.x+","+point.y+")");
		}*/
		for (Edge edge1 : edges) {

			int cp1 = compter(points, edge1.p);
			int cp2 = compter(points, edge1.q);
			
			if (cp1 == 1) {
				points.remove(edge1.p);
			}
			if (cp2 == 1) {
				points.remove(edge1.q);
			}
		}
		return points;	
	}
		
	
	private ArrayList<Point> eclatement_edges(ArrayList<Edge> edgesAlgo) {
		ArrayList<Point> pts = new ArrayList<>();
		for (Edge edge : edgesAlgo) {
			pts.add(edge.p);
			pts.add(edge.q);
		}
		return pts;
	}

	public Tree2D calculSteinerBudget(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
		ArrayList<Object> path_dist = calculShortestPaths(points, edgeThreshold);
		int[][] paths = (int[][]) path_dist.get(0);
		double[][] dist = (double[][]) path_dist.get(1);

		HashMap<Point, Integer> mapPointIndice = new HashMap<>();

		for (int i = 0; i < points.size(); i++) {
			mapPointIndice.put(points.get(i), i);
		}

		ArrayList<Edge> kruskhaled = kruskal(hitPoints, dist, mapPointIndice);

		ArrayList<Edge> best_edges = new ArrayList<>();
		for (Edge edge : kruskhaled) {
			Point p = edge.p;
			Point q = edge.q;

			ArrayList<Edge> edges_intermediaires = new ArrayList<>();

			int i = mapPointIndice.get(p);
			int j = mapPointIndice.get(q);

			if (paths[i][j] != j) {

				int k = paths[i][j];
				edges_intermediaires
						.add(new Edge(getKeyFromValue(mapPointIndice, i), getKeyFromValue(mapPointIndice, k)));
				while (k != j) {
					int k_moins = k;
					k = paths[k][j];
					edges_intermediaires.add(
							new Edge(getKeyFromValue(mapPointIndice, k_moins), getKeyFromValue(mapPointIndice, k)));

				}
				edges_intermediaires
						.add(new Edge(getKeyFromValue(mapPointIndice, k), getKeyFromValue(mapPointIndice, j)));
				best_edges.addAll(edges_intermediaires);

			} else {
				best_edges.add(edge);
			}

		}
		//System.out.println(best_edges.size());
	   //ArrayList<Edge> best_edges_filtres = filtrage_edges_extrems(filtrage_edges_extrems(filtrage_edges_extrems(filtrage_edges_extrems(best_edges))));
		ArrayList<Edge> best_edges_filtres = best_edges;
		int ds = distances_total_edges(best_edges_filtres);
		while(ds>2106){
		   best_edges_filtres = filtrage_edges_extrems(best_edges_filtres);
		   ds = distances_total_edges(best_edges_filtres);
		   System.out.println(" number of edges "+best_edges.size());
		   System.out.println(" distance "+ds);;
	   }
		
		ArrayList<Point> hitsPointsNew = filtrage_sommets_extrems(best_edges_filtres);
		ArrayList<Edge> kruskaled = kruskal(hitsPointsNew, dist, mapPointIndice);

		ArrayList<Edge> best_edges_new = new ArrayList<>();
		for (Edge edge : kruskaled) {
			Point p = edge.p;
			Point q = edge.q;

			ArrayList<Edge> edges_intermediaires = new ArrayList<>();

			int i = mapPointIndice.get(p);
			int j = mapPointIndice.get(q);

			if (paths[i][j] != j) {

				int k = paths[i][j];
				edges_intermediaires
						.add(new Edge(getKeyFromValue(mapPointIndice, i), getKeyFromValue(mapPointIndice, k)));
				while (k != j) {
					int k_moins = k;
					k = paths[k][j];
					edges_intermediaires.add(
							new Edge(getKeyFromValue(mapPointIndice, k_moins), getKeyFromValue(mapPointIndice, k)));

				}
				edges_intermediaires
						.add(new Edge(getKeyFromValue(mapPointIndice, k), getKeyFromValue(mapPointIndice, j)));
				best_edges_new.addAll(edges_intermediaires);

			} else {
				best_edges_new.add(edge);
			}

		}
		ArrayList<Edge> best_edges_filtres_new = best_edges_new;
		int ds1 = distances_total_edges(best_edges_filtres_new);
		while(ds1>1664){
		   best_edges_filtres_new = filtrage_edges_extrems(best_edges_filtres_new);
		   ds1 = distances_total_edges(best_edges_filtres_new);
		   System.out.println(" number of edges "+best_edges_new.size());
		   System.out.println(" distance "+ds1);;
	   }
	   //System.out.println(best_edges_filtres.size());
	   return edgesToTree(best_edges_filtres_new,hitPoints.get(0));
	  }

	private int distances_total_edges(ArrayList<Edge> best_edges) {
		int res = 0;
		for (Edge edge : best_edges) {
			res +=edge.distance();
		}
		return res;
	}

}

class Edge {
	protected Point p, q;

	protected Edge(Point p, Point q) {
		this.p = p;
		this.q = q;
	}

	protected double distance() {
		return p.distance(q);
	}
}

class NameTag {
	private ArrayList<Point> points;
	private int[] tag;

	@SuppressWarnings("unchecked")
	protected NameTag(ArrayList<Point> points) {
		this.points = (ArrayList<Point>) points.clone();
		tag = new int[points.size()];
		for (int i = 0; i < points.size(); i++)
			tag[i] = i;
	}

	protected void reTag(int j, int k) {
		for (int i = 0; i < tag.length; i++)
			if (tag[i] == j)
				tag[i] = k;
	}

	protected int tag(Point p) {
		for (int i = 0; i < points.size(); i++)
			if (p.equals(points.get(i)))
				return tag[i];
		return 0xBADC0DE;
	}
}
