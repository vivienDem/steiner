package algorithms;

import java.awt.Point;
import java.util.ArrayList;

public class Tree2D {
	private Point root;
	private ArrayList<Tree2D> subtrees;

	/**
	 * Constructeur d'un arbre
	 * 
	 * @param p     la racine
	 * @param trees la liste des enfants
	 */
	public Tree2D(Point p, ArrayList<Tree2D> trees) {
		this.root = p;
		this.subtrees = trees;
	}

	/**
	 * Accesseur sur la racine
	 * 
	 * @return la racine
	 */
	public Point getRoot() {
		return this.root;
	}

	/**
	 * Accesseur sur les enfants
	 * 
	 * @return les enfants
	 */
	public ArrayList<Tree2D> getSubTrees() {
		return this.subtrees;
	}

	/**
	 * Calcule la distance entre la racine et tous les enfants
	 * 
	 * @return la distance entre la racine et tous les enfants
	 */
	public double distanceRootToSubTrees() {
		double d = 0;
		for (int i = 0; i < this.subtrees.size(); i++)
			d += subtrees.get(i).getRoot().distance(root);
		return d;
	}
}
