package algorithms;

import java.awt.Point;
import java.util.ArrayList;

public class NameTag {
	private ArrayList<Point> points;
	private int[] tag;

	/**
	 * Constructeur d'etiquettes
	 * 
	 * @param points la liste des points a etiquetter
	 */
	protected NameTag(ArrayList<Point> points) {
		this.points = (ArrayList<Point>) points.clone();
		tag = new int[points.size()];
		for (int i = 0; i < points.size(); i++) {
			tag[i] = i;
		}
	}

	/**
	 * Renommage d'etiquettes
	 * 
	 * @param j l'etiquette a renommer
	 * @param k la nouvelle valeur
	 */
	protected void reTag(int j, int k) {
		for (int i = 0; i < tag.length; i++) {
			if (tag[i] == j)
				tag[i] = k;
		}
	}

	/**
	 * Accesseur sur le tag du point p
	 * 
	 * @param p le point dont on cherche le tag
	 * @return le tag du point
	 */
	protected int tag(Point p) {
		for (int i = 0; i < points.size(); i++) {
			if (p.equals(points.get(i)))
				return tag[i];
		}

		return 0xBADC0DE;
	}
}
