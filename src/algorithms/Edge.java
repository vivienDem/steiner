package algorithms;

import java.awt.Point;
import java.util.Comparator;

public class Edge implements Comparable<Edge> {
	protected Point p, q;

	/**
	 * Constructeur d'une arete
	 * 
	 * @param p le premier point
	 * @param q le deuxieme point
	 */
	protected Edge(Point p, Point q) {
		this.p = p;
		this.q = q;
	}

	/**
	 * Verifie si un point appartient a l'arete
	 * 
	 * @param x le point
	 * @return true si x appartient a l'arete, false sinon
	 */
	public boolean contains(Point x) {
		return p == x || q == x;
	}

	/**
	 * Accesseur du premier point de l'arete
	 * 
	 * @return le premier point
	 */
	public Point getP1() {
		return p;
	}

	/**
	 * Accesseur du deuxieme point de l'arete
	 * 
	 * @return le deuxieme point
	 */
	public Point getP2() {
		return q;
	}

	/**
	 * Accesseur du poids de l'arete
	 * 
	 * @return le poids
	 */
	public double getWeight() {
		return p.distance(q);
	}

	public double getAngleWith(Edge e2) {
		Point A = CommonPointWith(e2);
		Point B = UncommonPointWith(e2);
		Point C = e2.UncommonPointWith(this);
		if (A == null || B == null || C == null)
			return 0;
		return angle(A, B, A, C);
	}

	private double angle(Point p, Point q, Point s, Point t) {
		if (p.equals(q) || s.equals(t))
			return Double.POSITIVE_INFINITY;
		double cosTheta = ((q.x - p.x) * (t.x - s.x) + (q.y - p.y) * (t.y - s.y))
				/ (double) (p.distance(q) * s.distance(t));
		return Math.acos(cosTheta);
	}

	private Point CommonPointWith(Edge e2) {
		if ((p == e2.getP1() && q != e2.getP2()) || (q != e2.getP1() && p == e2.getP2()))
			return p;
		if ((p != e2.getP1() && q == e2.getP2()) || (q == e2.getP1() && p != e2.getP2()))
			return q;
		return null;
	}

	private Point UncommonPointWith(Edge e2) {
		if ((p == e2.getP1() && q != e2.getP2()) || (q != e2.getP1() && p == e2.getP2()))
			return q;
		if ((p != e2.getP1() && q == e2.getP2()) || (q == e2.getP1() && p != e2.getP2()))
			return p;
		return null;
	}

	@Override
	public String toString() {
		return p.toString() + " " + q.toString() + " weight : " + getWeight();

	}

	@Override
	public int compareTo(Edge o) {
		return Double.compare(getWeight(), o.getWeight());
	}

}
