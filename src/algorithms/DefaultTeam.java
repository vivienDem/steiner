package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class DefaultTeam {

	/**
	 * Steiner sans budget
	 * 
	 * @param points        la liste des points
	 * @param edgeThreshold le seuil
	 * @param hitPoints     la liste des points a relier
	 * @return un arbre couvrant les hitPoints
	 */
	public Tree2D calculSteiner(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
		int[][] matrice_acces = calculShortestPaths(points, edgeThreshold);
		ArrayList<Edge> edges = calculSteinerIntermediaire(points, matrice_acces, hitPoints);

		ArrayList<Edge> res = (ArrayList<Edge>) edges.clone();
		Point A, B, C, I;
		Point barycenter, fermat;
		double score, scoreBarycenter, scoreFermat, bestScore, oldScore, newScore = 0;
		oldScore = getDistance(edges);
		int cpt = 0, size;
		while (newScore < oldScore) {
			cpt++;
			size = edges.size();
			for (int i = 0; i < size; i++) {
				for (int j = i + 1; j < size; j++) {
					Edge e1 = edges.get(i);
					Edge e2 = edges.get(j);

					if (e1.getP1().equals(e2.getP1())) {
						A = e1.getP2();
						B = e1.getP1();
						C = e2.getP2();
					} else if (e1.getP1().equals(e2.getP2())) {
						A = e1.getP2();
						B = e1.getP1();
						C = e2.getP1();
					} else if (e2.getP1().equals(e1.getP2())) {
						A = e1.getP1();
						B = e1.getP2();
						C = e2.getP2();
					} else if (e1.getP2().equals(e2.getP2())) {
						A = e1.getP1();
						B = e1.getP2();
						C = e2.getP1();
					} else
						continue;

					barycenter = calculBarycenter(A, B, C);
					fermat = calculFermat(A, B, C);
					score = A.distance(B) + B.distance(C);
					scoreBarycenter = A.distance(barycenter) + B.distance(barycenter) + C.distance(barycenter);
					scoreFermat = A.distance(fermat) + B.distance(fermat) + C.distance(fermat);

					if (scoreFermat < scoreBarycenter) {
						bestScore = scoreFermat;
						I = fermat;
					} else {
						bestScore = scoreBarycenter;
						I = barycenter;
					}

					I = closestPoint(I, points);

					if (bestScore < score) {
						hitPoints.add(I);
						edges.add(new Edge(A, I));
						edges.add(new Edge(B, I));
						edges.add(new Edge(C, I));
						edges.remove(e1);
						edges.remove(e2);
						i--;
						break;
					}
				}
			}
			edges = calculSteinerIntermediaire(points, matrice_acces, hitPoints);
			newScore = getDistance(edges);
			res = edges;
			oldScore = newScore;
		}
		return edgesToTree(res, res.get(0).getP1());
	}

	/**
	 * Steiner avec un budget de 1664
	 * 
	 * @param points        la liste des points
	 * @param edgeThreshold le seuil
	 * @param hitPoints     la liste des points a relier
	 * @return un arbre couvrant les hitPoints
	 */
	public Tree2D calculSteinerBudget(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
		int[] indices = new int[hitPoints.size()];
		int i = 0, j;
		for (Point p : hitPoints) {
			j = 0;
			for (Point q : points) {
				if (p == q)
					indices[i] = j;
				j++;
			}
			i++;
		}
		ArrayList<Edge> edges = new ArrayList<Edge>();
		Tree2D K = calculPrim(hitPoints);
		int[][] matrice_acces = calculShortestPaths(points, edgeThreshold);
		ArrayList<Tree2D> trees = new ArrayList<Tree2D>();
		trees.add(K);
		while (!trees.isEmpty()) {
			Tree2D tmp = trees.remove(0);
			Point root = tmp.getRoot();
			ArrayList<Tree2D> children = tmp.getSubTrees();
			int root_index = points.indexOf(root);
			for (Tree2D e : children) {
				trees.add(e);
				int index_child = points.indexOf(e.getRoot());
				int index = root_index;
				while (index != index_child) {
					Point a = points.get(index);
					Point b = points.get(matrice_acces[index][index_child]);
					Edge edge = new Edge(a, b);
					edges.add(edge);
					index = matrice_acces[index][index_child];
				}
			}
		}
		double distance = getDistance(edges);
		ArrayList<Edge> bordures = getBorders(edges);
		while (distance > 1664) {
			bordures = getBorders(edges);
			Edge e = removeMaxEdge(bordures);
			edges.remove(e);
			distance -= e.getWeight();
		}
		return edgesToTree(edges, K.getRoot());
	}

	/**
	 * Enleve l'arete la plus lourde
	 * 
	 * @param bordures la liste des aretes au bord
	 * @return l'arete qui a ete enlevee
	 */
	private Edge removeMaxEdge(ArrayList<Edge> bordures) {
		Edge enlever = null;
		double weight = 0;
		for (Edge e : bordures) {
			if (e.getWeight() > weight) {
				weight = e.getWeight();
				enlever = e;
			}
		}
		bordures.remove(enlever);
		return enlever;
	}

	/**
	 * Calcule la distance totale des aretes
	 * 
	 * @param edges la liste des aretes
	 * @return la distance totale des aretes
	 */
	public double getDistance(ArrayList<Edge> edges) {
		double res = 0;
		for (Edge e : edges) {
			res += e.getWeight();
		}
		return res;
	}

	/**
	 * Construit la liste des aretes au bord
	 * 
	 * @param edges la liste des aretes
	 * @return la liste des aretes au bord
	 */
	public ArrayList<Edge> getBorders(ArrayList<Edge> edges) {
		ArrayList<Edge> res = new ArrayList<Edge>();
		for (Edge edge : edges) {
			if (border(edge.getP1(), edge.getP2(), edges))
				res.add(edge);
		}
		return res;
	}

	/**
	 * Calcule si une arete relie un noeud interne a une feuille
	 * 
	 * @param p     un point de l'arete
	 * @param q     l'autre point de l'arete
	 * @param edges la liste des aretes
	 * @return true si l'arete est au bord, false sinon
	 */
	public boolean border(Point p, Point q, ArrayList<Edge> edges) {
		int nbP = 0, nbQ = 0;
		for (Edge e : edges) {
			if (e.contains(p))
				nbP++;
			if (e.contains(q))
				nbQ++;
			if (nbP > 1 && nbQ > 1)
				return false;
		}
		return true;
	}

	/**
	 * Transforme une liste d'aretes en arbre
	 * 
	 * @param edges la liste d'aretes
	 * @param root  la racine
	 * @return un arbre de racine root contenant les aretes
	 */
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

	/**
	 * Calcule les chemins les plus courts entre chaque points
	 * 
	 * @param points        la liste des points
	 * @param edgeThreshold le seuil a partir du quel deux points ne sont pas
	 *                      adjacents
	 * @return une matrice d'acces
	 */
	public int[][] calculShortestPaths(ArrayList<Point> points, int edgeThreshold) {
		int[][] paths = new int[points.size()][points.size()];
		double[][] distances = new double[points.size()][points.size()];
		double distance;
		int i = 0, j;
		for (i = 0; i < paths.length; i++) {
			for (j = 0; j < paths.length; j++) {
				if (i == j) {
					distances[i][i] = 0;
					continue;
				}
				distance = points.get(i).distance(points.get(j));
				if (distance <= edgeThreshold) {
					distances[i][j] = distance;
				} else {
					distances[i][j] = Double.POSITIVE_INFINITY;
				}
				paths[i][j] = j;
			}
		}

		for (int k = 0; k < points.size(); k++) {
			for (i = 0; i < points.size(); i++) {
				for (j = 0; j < points.size(); j++) {
					if (i == j)
						continue;
					if (distances[i][k] + distances[k][j] < distances[i][j]) {

						distances[i][j] = distances[i][k] + distances[k][j];
						paths[i][j] = paths[i][k];
					}
				}
			}
		}
		return paths;
	}

	/**
	 * Trouve le point le plus proche de point dans la liste points
	 * 
	 * @param point  le point dont on cherche son voisin le plus proche
	 * @param points la liste des points
	 * @return le point appartenant a points le plus proche de point
	 */
	private Point closestPoint(Point point, ArrayList<Point> points) {
		double distance = Double.POSITIVE_INFINITY;
		Point res = null;
		for (Point p : points) {
			if (point.distance(p) < distance) {
				distance = point.distance(p);
				res = p;
			}
		}
		return res;
	}

	/**
	 * Calcule le point de Fermat a partir des trois points A, B et C
	 * 
	 * @param A le premier point
	 * @param B le second point
	 * @param C le troisieme point
	 * @return le point de Fermat
	 */
	private Point calculFermat(Point A, Point B, Point C) {
		Edge AB = new Edge(A, B);
		Edge BC = new Edge(B, C);
		Edge AC = new Edge(A, C);

		if (Math.abs(AB.getAngleWith(BC)) >= 120)
			return B;
		else if (Math.abs(AB.getAngleWith(AC)) >= 120)
			return A;
		else if (Math.abs(BC.getAngleWith(AC)) >= 120)
			return C;

		Point E = null, F = null;
		int vectABx = B.x - A.x;
		int vectABy = B.y - A.y;
		int vectACx = C.x - A.x;
		int vectACy = C.y - A.y;

		double lcos = Math.cos(Math.PI / 3);
		double rcos = Math.cos(-Math.PI / 3);
		double lsin = Math.sin(Math.PI / 3);
		double rsin = Math.sin(-Math.PI / 3);
		if ((vectABx * vectACy - vectABy * vectACx) > 0) {
			E = new Point((int) Math.round(A.x + vectACx * lcos - vectACy * lsin),
					(int) Math.round(A.y + vectACy * lcos + vectACx * lsin));
			F = new Point((int) Math.round(A.x + vectABx * rcos - vectABy * rsin),
					(int) Math.round(A.y + vectABy * rcos + vectABx * rsin));
		} else {
			E = new Point((int) Math.round(A.x + vectACx * rcos - vectACy * rsin),
					(int) Math.round(A.y + vectACy * rcos + vectACx * rsin));
			F = new Point((int) Math.round(A.x + vectABx * lcos - vectABy * lsin),
					(int) Math.round(A.y + vectABy * lcos + vectABx * lsin));
		}

		double A1 = C.y - F.y;
		double B1 = F.x - C.x;
		double C1 = A1 * F.x + B1 * F.y;

		double A2 = B.y - E.y;
		double B2 = E.x - B.x;
		double C2 = A2 * E.x + B2 * E.y;

		double det = A1 * B2 - A2 * B1;
		double x, y;
		if (det == 0) {
			return A;
		} else {
			x = (B2 * C1 - B1 * C2) / det;
			y = (A1 * C2 - A2 * C1) / det;
		}
		return new Point((int) x, (int) y);
	}

	/**
	 * Calcule le barycentre a partir des trois points a, b et c
	 * 
	 * @param a le premier point
	 * @param b le second point
	 * @param c le troisieme point
	 * @return le barycentre
	 */
	private Point calculBarycenter(Point a, Point b, Point c) {
		return new Point((int) (a.getX() + b.getX() + c.getX()) / 3, (int) (a.getY() + b.getY() + c.getY()) / 3);
	}

	/**
	 * Construit l'arbre de Steiner en utilisant Kruskal et une matrice d'acces
	 * 
	 * @param points        la liste des points
	 * @param matrice_acces une matrice contenant les chemins pour passer d'un point
	 *                      a l'autre
	 * @param hitPoints     les points a relier
	 * @return un arbre couvrant les hitPoints
	 */
	public ArrayList<Edge> calculSteinerIntermediaire(ArrayList<Point> points, int[][] matrice_acces,
			ArrayList<Point> hitPoints) {

		ArrayList<Edge> spanningTreeEdges = calculKruskal(hitPoints);
		ArrayList<Point> path = new ArrayList<>();
		for (Edge e : spanningTreeEdges) {
			int i = points.indexOf(e.getP1());
			int j = points.indexOf(e.getP2());
			ArrayList<Integer> pointsIJ = getPath(i, j, matrice_acces);
			for (Integer k : pointsIJ) {
				path.add(points.get(k));
			}
		}
		return calculKruskal(path);

	}

	/**
	 * Trouve le cout minimum des sommets non visites
	 * 
	 * @param cost    un tableau des couts des sommets
	 * @param reached un tableau indiquant si un sommet a ete visite
	 * @return l'indice du sommet avec le cout minimum
	 */
	public int getMinCost(double cost[], boolean reached[]) {
		double min = Double.POSITIVE_INFINITY;
		int min_index = 0;
		for (int i = 0; i < cost.length; i++) {
			if (!reached[i] && cost[i] < min) {
				min = cost[i];
				min_index = i;
			}
		}
		return min_index;
	}

	/**
	 * Algorithme de Kruskal
	 * 
	 * @param points la liste des points
	 * @return l'arbre couvrant minimum
	 */
	public ArrayList<Edge> calculKruskal(ArrayList<Point> points) {
		List<Edge> edges = new ArrayList<Edge>();
		for (Point p : points) {
			for (Point q : points) {
				if (p.equals(q) || contains(edges, p, q))
					continue;
				edges.add(new Edge(p, q));
			}
		}

		Collections.sort(edges);

		ArrayList<Edge> kruskal = new ArrayList<Edge>();
		Edge current;
		NameTag forest = new NameTag(points);
		while (edges.size() != 0) {
			current = edges.remove(0);
			if (forest.tag(current.p) != forest.tag(current.q)) {
				kruskal.add(current);
				forest.reTag(forest.tag(current.p), forest.tag(current.q));
			}
		}

		return kruskal;
	}

	/**
	 * Algorithme de Prim
	 * 
	 * @param points la liste des points
	 * @return l'arbre couvrant minimum
	 */
	public Tree2D calculPrim(ArrayList<Point> points) {
		double[] cost = new double[points.size()];
		int[] pred = new int[points.size()];
		boolean[] reached = new boolean[points.size()];
		ArrayList<Edge> edges = new ArrayList<Edge>();

		for (int i = 0; i < cost.length; i++) {
			cost[i] = Double.POSITIVE_INFINITY;
			reached[i] = false;
		}
		cost[0] = 0;
		pred[0] = -1;

		for (int i = 0; i < points.size() - 1; i++) {
			int j = getMinCost(cost, reached);
			reached[j] = true;

			for (int k = 0; k < points.size(); k++) {
				double distance = points.get(j).distance(points.get(k));
				if (distance != 0 && !reached[k] && distance < cost[k]) {
					pred[k] = j;
					cost[k] = distance;
				}
			}
		}

		for (int i = 1; i < points.size(); i++) {
			Point p = points.get(i);
			Point q = points.get(pred[i]);
			edges.add(new Edge(p, q));
		}
		return edgesToTree(edges, edges.get(0).getP1());
	}

	/**
	 * Verifie si edges contient une arete composee de p et q
	 * 
	 * @param edges la liste des aretes
	 * @param p     le premier point
	 * @param q     le deuxieme point
	 * @return true si le predicat est verifie, false sinon
	 */
	private boolean contains(List<Edge> edges, Point p, Point q) {

		for (Edge e : edges) {
			if (e.p.equals(p) && e.q.equals(q) || e.p.equals(q) && e.q.equals(p))
				return true;
		}
		return false;
	}

	/**
	 * Trouve le chemin entre les points i et j
	 * 
	 * @param i     le premier point
	 * @param j     le deuxieme point
	 * @param paths la matrice d'acces
	 * @return la liste des indices des points qu'il faut parcourir pour aller de i
	 *         vers j
	 */
	public static ArrayList<Integer> getPath(int i, int j, int[][] paths) {
		ArrayList<Integer> path = new ArrayList<>();
		path.add(i);
		while (paths[i][j] != j) {
			path.add(paths[i][j]);
			i = paths[i][j];
		}
		path.add(j);
		return path;
	}
}
