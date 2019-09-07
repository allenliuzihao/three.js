/**
 *	@author zz85 / http://twitter.com/blurspline / http://www.lab4games.net/zz85/blog
 *	@author centerionware / http://www.centerionware.com
 *  @author Zihao Liu / catmull-clark subdivision and texture mapping onto the surface
 *
 *	Subdivision Geometry Modifier
 *		using Catmull-Clark and Loop Subdivision Scheme
 *
 *	References:
 *  	Catmull-Clark Subdivision:
 * 			http://users.cms.caltech.edu/~cs175/cs175-02/resources/CC.pdf
 * 			https://en.wikipedia.org/wiki/Catmull%E2%80%93Clark_subdivision_surface
 *		Loop Subdivision
 *			http://graphics.stanford.edu/~mdfisher/subdivision.html
 *			http://www.holmes3d.net/graphics/subdivision/
 *			http://www.cs.rutgers.edu/~decarlo/readings/subdiv-sg00c.pdf
 *
 *	Known Issues:
 *		- currently doesn't handle "Sharp Edges"
 */

import { Face3 } from '../../../src/core/Face3.js';
import { Quad } from '../../../src/core/Quad.js';
import { Geometry } from '../../../src/core/Geometry.js';
import { Vector2 } from '../../../src/math/Vector2.js';
import { Vector3 } from '../../../src/math/Vector3.js';

const CATMULL_CLARK_SCHEME = "catmull_clark";
const LOOP_SCHEME = "loop";

var SubdivisionModifier = function ( subdivisions, scheme ) {

	this.subdivisions = ( subdivisions === undefined ) ? 1 : subdivisions;

	this.scheme = scheme;

};

// Applies the "modify" pattern
/**
 * @param {Geometry} geometry
 */
SubdivisionModifier.prototype.modify = function ( geometry ) {

	if ( geometry.isBufferGeometry ) {

		geometry = new Geometry().fromBufferGeometry( geometry );

	} else {

		geometry = geometry.clone();

	}

	geometry.mergeVertices();

	var repeats = this.subdivisions;

	while ( repeats -- > 0 ) {

		this.smooth( geometry );

	}

	if(this.scheme === CATMULL_CLARK_SCHEME){
		geometry.toTriangleMesh();
	}

	geometry.computeFaceNormals();
	geometry.computeVertexNormals();

	return geometry;

};

( function () {

	/////////////////////////////
	// Performs one iteration of Subdivision
	SubdivisionModifier.prototype.smooth = function ( geometry ) {
		if(this.scheme === LOOP_SCHEME) {
			Loop_Subdivision (geometry);
		} else if (this.scheme === CATMULL_CLARK_SCHEME){
			Catmull_Clark_Subdivision (geometry);
		}
	};

	/* Catmull-Clark subdivision scheme */
	const EDGE_SPLIT = "_";

	/**
	 * Given an array of vectors, compute its average
	 * @param {Vector3[]} vectors 
	 */
	function avg(vectors){
		if(vectors.length === 0){
			return new Vector3();
		}

		let result = new Vector3();
		for (const vector of vectors) {
			result.add(vector);
		}
		return result.multiplyScalar(1/vectors.length);
	}

	/**
	 * generate edge key given vertex indices of the edge.
	 * @param {Number} vertexAIdx 
	 * @param {Number} vertexBIdx 
	 */
	function generateEdgeKey (vertexAIdx, vertexBIdx) {
		let vertexIndexA = Math.min( vertexAIdx, vertexBIdx );
		let vertexIndexB = Math.max( vertexAIdx, vertexBIdx );
		return vertexIndexA + EDGE_SPLIT + vertexIndexB;
	}

	/**
	 * generate "vertex1-vertex2" to faces set per edge
	 * 
	 * @param {Number} vertexAIdx 
	 * @param {Number} vertexBIdx 
	 * @param {Face3|Quad} face 
	 * @param {Map<String, Set<Face3|Quad>>} edgeToFaces 
	 */
	function generateEachEdge (vertexAIdx, vertexBIdx, face, edgeToFaces){
		let key = generateEdgeKey(vertexAIdx, vertexBIdx);

		if(!(key in edgeToFaces)){
			edgeToFaces.set(key, new Set(face));
		} else {
			edgeToFaces.get(key).add(face);
		}
	}

	/**
	 * generate "vertex1-vertex2" to faces set per face
	 * 
	 * @param {Face3[]|Quad[]} faces 
	 * @param {Map<String, Set<Face3|Quad>>} edgeToFaces 
	 */
	function generateEdges(faces, edgeToFaces) {
		let face;

		for (let i = 0, il = faces.length; i < il; i++ ) {
			face = faces[i];
			generateEachEdge(face.a, face.b, face, edgeToFaces);
			generateEachEdge(face.b, face.c, face, edgeToFaces);
			generateEachEdge(face.c, face.a, face, edgeToFaces);
		}
	}

	/**
	 * add a face information to a vertex.
	 * @param {Number} vertexIdx 
	 * @param {Face3|Quad} face 
	 * @param {Map<Number, Set<Face3|Quad>>} vertexToFaces 
	 */
	function addFacePerVertex (vertexIdx, face, vertexToFaces){
		let faceSet;
		if(vertexIdx in vertexToFaces){
			faceSet = vertexToFaces.get(vertexIdx);
		} else {
			faceSet = new Set();
			vertexToFaces.set(vertexIdx, faceSet);
		}
		faceSet.add(face);
	}

	/**
	 * add edge information to a vertex.
	 * @param {Number} vertexAIdx 
	 * @param {Number} vertexBIdx 
	 * @param {Map<Number, Set<String>>} vertexToEdges
	 */
	function addEdgePerVertex (vertexAIdx, vertexBIdx, vertexToEdges){
		let edgeSet;
		if(vertexAIdx in vertexToEdges){
			edgeSet = vertexToEdges.get(vertexAIdx);
		} else {
			edgeSet = new Set();
			vertexToEdges.set(vertexAIdx, edgeSet);
		}
		edgeSet.add(generateEdgeKey(vertexAIdx, vertexBIdx));

		if(vertexBIdx in vertexToEdges){
			edgeSet = vertexToEdges.get(vertexBIdx);
		} else {
			edgeSet = new Set();
			vertexToEdges.set(vertexBIdx, edgeSet);
		}
		edgeSet.add(generateEdgeKey(vertexAIdx, vertexBIdx));
	}

	/**
	 * generate mapping of vertex -> {set {faces}, set {edges}}
	 * basically all the incident faces and edges for that vertex 
	 * 
	 * @param {Face3[]|Quad[]} faces 
	 * @param {Map<Number, Set<Face3|Quad>>} vertexToFaces
	 * @param {Map<Number, Set<String>>} vertexToEdges
	 */
	function generateVertexInfo(faces, vertexToFaces, vertexToEdges){
		for (const face of faces) {
			if(face instanceof Face3){
				addFacePerVertex(face.a, vertexToFaces);
				addFacePerVertex(face.b, vertexToFaces);
				addFacePerVertex(face.c, vertexToFaces);

				addEdgePerVertex(face.a, face.b, vertexToEdges);
				addEdgePerVertex(face.b, face.c, vertexToEdges);
				addEdgePerVertex(face.c, face.a, vertexToEdges);
			} else if (face instanceof Quad){
				addFacePerVertex(face.a, vertexToFaces);
				addFacePerVertex(face.b, vertexToFaces);
				addFacePerVertex(face.c, vertexToFaces);
				addFacePerVertex(face.d, vertexToFaces);

				addEdgePerVertex(face.a, face.b, vertexToEdges);
				addEdgePerVertex(face.b, face.c, vertexToEdges);
				addEdgePerVertex(face.c, face.d, vertexToEdges);
				addEdgePerVertex(face.d, face.a, vertexToEdges);
			}
		}
	}

	/**
	 * calculate face point per face
	 * 
	 * @param {Face3[]|Quad[]} faces 
	 * @param {Vector3[]} vertices 
	 * @param {Map<Face3|Quad, Vector3>} faceToFacePoints 
	 */
	function generateFacePoints (faces, vertices, faceToFacePoints){
		let face, vertexA, vertexB, vertexC;

		for (let index = 0, length = faces.length; index < length; index++) {
			face = faces[index];

			if (face instanceof Face3){
				vertexA = vertices[face.a];
				vertexB = vertices[face.b];
				vertexC = vertices[face.c];
				faceToFacePoints.set(face, 
					new Vector3()
						.addVectors(vertexA, vertexB)
						.add(vertexC)
						.multiplyScalar(1/3));
			} else if (face instanceof Quad){
				vertexA = vertices[face.a];
				vertexB = vertices[face.b];
				vertexC = vertices[face.c];
				let vertexD = vertices[face.d];
				faceToFacePoints.set(face, 
					new Vector3()
						.addVectors(vertexA, vertexB)
						.addVectors(vertexC, vertexD)
						.multiplyScalar(1/4));
			}
		}
	}

	/**
	 * for each edge, generate an edge point.
	 * 
	 * @param {Vector3[]} vertices 
	 * @param {Map<String, Set<Face3|Quad>>} edgeToFaces
	 * @param {Map<Face3|Quad, Vector3>} faceToFacePoints
	 * @param {Map<String, Vector3>} edgeToEdgePoints
	 */
	function generateEdgePoints (vertices, edgeToFaces, faceToFacePoints, edgeToEdgePoints) {
		let edgeSplits, facesArray;
		let vertexA, vertexB, facePointA, facePointB;
		let edgePoint;
	
		for (let [edge, faces] of edgeToFaces) {
			edgeSplits = edge.split(EDGE_SPLIT);
			vertexA = vertices[Number(edgeSplits[0])];
			vertexB = vertices[Number(edgeSplits[1])];

			if(faces.size !== 2){
				console.error("SubdivisionModifier generateEdgePoints error: each edge should only have two incident faces.");
			}

			facesArray = Array.from(faces);
			facePointA = faceToFacePoints.get(facesArray[0]);
			facePointB = faceToFacePoints.get(facesArray[1]);

			edgePoint = 
				new Vector3()
				.addVectors(vertexA, vertexB)
				.addVectors(facePointA, facePointB)
				.multiplyScalar(1/4);
			
			edgeToEdgePoints.set(edge, edgePoint);
		}
	}

	/**
	 * 
	 * @param {Vector3[]} vertices 
	 * @param {Map<Number, Set<Face3|Quad>>} vertexToFaces 
	 * @param {Map<Number, Set<String>>} vertexToEdges 
	 * @param {Map<Face3|Quad, Vector3>} faceToFacePoint
	 * @param {Vector3[]} nVertices 
	 */
	function generateNewVertexPosition (vertices, vertexToFaces, vertexToEdges, faceToFacePoint, nVertices) {
		if(nVertices.length !== vertices.length){
			return;
		}

		let vertex, avgFacePoints, avgEdgeMidPoints, n;
		let edgeSplits, vertexA, vertexB;
		let faces, edges, facePoints, edgeMidPoints;
		for (let index = 0; index < vertices.length; index++) {
			vertex = vertices[index];
			faces = vertexToFaces.get(index);
			edges = vertexToEdges.get(index);
			facePoints = [];
			edgeMidPoints = [];

			if(faces.size == 0 || edges.size == 0){
				console.error("SubdivisionModifier generateNewVertexPosition: faces and edges for vertex should not be empty.", this);
			}

			for (const face of faces) {
				if(face in faceToFacePoint){
					facePoints.push(faceToFacePoint.get(face));
				} else {
					console.error("SubdivisionModifier generateNewVertexPosition: face should have a face point.", this);
				}
			}

			for (const edge of edges) {
				edgeSplits = edge.split(EDGE_SPLIT);
				vertexA = vertices[Number(edgeSplits[0])];
				vertexB = vertices[Number(edgeSplits[1])];
				edgeMidPoints.push(new Vector3().addVectors(vertexA, vertexB).multiplyScalar(2.0));
			}

			if(facePoints.length === 0 || edgeMidPoints.length === 0){
				console.error("SubdivisionModifier generateNewVertexPosition: number of face points and edge midpoints should not be empty.", this);
			}

			if(facePoints.length !== edgeMidPoints.length) {
				console.error("SubdivisionModifier generateNewVertexPosition: number of face points should match the number of edge midpoints", this);
			}

			n = facePoints.length;
			avgFacePoints = avg(facePoints);
			avgEdgeMidPoints = avg(edgeMidPoints);
			nVertices[index] = 
				new Vector3()
				.add(avgFacePoints)
				.addScaledVector(avgEdgeMidPoints, 2)
				.addScaledVector(vertex, n - 3)
				.multiplyScalar(1/n);
		}
	}

	/**
	 * generate quad faces and new vertices from face points and edge points.
	 * 
	 * @param {Vector3[]} vertices
	 * @param {Map<Face3|Quad, Vector3>} faceToFacePoint
	 * @param {Map<String, Vector3>} edgeToEdgePoint
	 * @returns object containing vertices and quads of smoothed mesh
	 */
	function generateSmoothedMesh (vertices, faceToFacePoint, edgeToEdgePoint) {
		let facePoints = new Array(), edgePoints = new Array();
		let quads = new Array();

		let edgePointIndex = 0, facePointIndex = 0;

		for (let [face, facePoint] of faceToFacePoint) {
			if (face instanceof Face3){
				facePoints.push(facePoint);
				edgePoints.push(edgeToEdgePoint.get(generateEdgeKey(face.a, face.b)));
				edgePoints.push(edgeToEdgePoint.get(generateEdgeKey(face.b, face.c)));
				edgePoints.push(edgeToEdgePoint.get(generateEdgeKey(face.c, face.a)));

				quads.push(
					new Quad(
						edgePointIndex + vertices.length + faceToFacePoint.size, 
						face.b, 
						edgePointIndex+ 1 + vertices.length + faceToFacePoint.size, 
						facePointIndex + vertices.length,
						undefined,
						undefined,
						face.materialIndex));
				quads.push(
					new Quad(
						facePointIndex + vertices.length,
						edgePointIndex + 1 + vertices.length + faceToFacePoint.size,
						face.c, 
						edgePointIndex + 2 + vertices.length + faceToFacePoint.size,
						undefined,
						undefined,
						face.materialIndex));
				quads.push(
					new Quad(
						facePointIndex + vertices.length,
						edgePointIndex + vertices.length + faceToFacePoint.size, 
						face.a, 
						edgePointIndex + 2 + vertices.length + faceToFacePoint.size,
						undefined,
						undefined,
						face.materialIndex));

				facePointIndex++;
				edgePointIndex += 3;
			} else if (face instanceof Quad){
				facePoints.push(facePoint);
				edgePoints.push(edgeToEdgePoint.get(generateEdgeKey(face.a, face.b)));
				edgePoints.push(edgeToEdgePoint.get(generateEdgeKey(face.b, face.c)));
				edgePoints.push(edgeToEdgePoint.get(generateEdgeKey(face.c, face.d)));
				edgePoints.push(edgeToEdgePoint.get(generateEdgeKey(face.d, face.a)));

				quads.push(
					new Quad(
						face.a,
						edgePointIndex + vertices.length + faceToFacePoint.size,
						facePointIndex + vertices.length, 
						edgePointIndex+ 3 + vertices.length + faceToFacePoint.size,
						undefined,
						undefined,
						face.materialIndex));
				quads.push(
					new Quad(
						edgePointIndex + vertices.length + faceToFacePoint.size,
						face.b, 
						edgePointIndex + 1 + vertices.length + faceToFacePoint.size,
						facePointIndex + vertices.length,
						undefined,
						undefined,
						face.materialIndex));
				quads.push(
					new Quad(
						facePointIndex + vertices.length,
						edgePointIndex + 1 + vertices.length + faceToFacePoint.size, 
						face.c, 
						edgePointIndex + 2 + vertices.length + faceToFacePoint.size,
						undefined,
						undefined,
						face.materialIndex));
				quads.push(
					new Quad(
						edgePointIndex + 3 + vertices.length + faceToFacePoint.size,
						facePointIndex + vertices.length, 
						edgePointIndex + 2 + vertices.length + faceToFacePoint.size,
						face.d,
						undefined,
						undefined,
						face.materialIndex));

				facePointIndex++;
				edgePointIndex += 4;
			}
		}

		return { "vertices": vertices.concat(facePoints).concat(edgePoints),"quads": quads };
	}

	/**
 	 * @param {Geometry} geometry
 	 */
	function Catmull_Clark_Subdivision (geometry){
		if(geometry === undefined || geometry.vertices === undefined
			|| geometry.faces == undefined || geometry.vertices.length === 0
			|| geometry.faces.length === 0){
				console.error("subdivision input geometry is not valid.", geometry);
				return;
			}

		let vertices = geometry.vertices, faces = geometry.faces;
		let edgeToFaces = new Map()
			,faceToFacePoint = new Map()
			,edgeToEdgePoint = new Map()
			,vertexToFaces = new Map()
			,vertexToEdges = new Map()
			,nVertices = new Array(vertices.length);

		// generate edges information; each edge information is "vertex1-vertex2" -> set {faces}
		generateEdges(faces, edgeToFaces);

		// generate vertex -> {set {faces} , set {edges}}
		generateVertexInfo(faces, edgeToFaces, vertexToFaces, vertexToEdges);

		// calculate the face points for all faces
		generateFacePoints(faces, vertices, faceToFacePoint);

		// calculate edge points
		generateEdgePoints(vertices, edgeToFaces, faceToFacePoint, edgeToEdgePoint);

		// calculate new vertex positions
		generateNewVertexPosition(vertices, vertexToFaces, vertexToEdges, faceToFacePoint, nVertices);

		// generate smoothed mesh
		let smoothed = generateSmoothedMesh(nVertices, faceToFacePoint, edgeToEdgePoint);

		// set the new geometry
		geometry.vertices = smoothed["vertices"];
		geometry.faces = smoothed["quads"];

		// TODO: vertex color and uv interpolation
		
	}

	/* Loop subdivision */
	
	// Some constants
	var ABC = [ 'a', 'b', 'c' ];
	
	function getEdge( a, b, map ) {

		var vertexIndexA = Math.min( a, b );
		var vertexIndexB = Math.max( a, b );

		var key = vertexIndexA + "_" + vertexIndexB;

		return map[ key ];

	}

	function processEdge( a, b, vertices, map, face, metaVertices ) {

		var vertexIndexA = Math.min( a, b );
		var vertexIndexB = Math.max( a, b );

		var key = vertexIndexA + "_" + vertexIndexB;

		var edge;

		if ( key in map ) {

			edge = map[ key ];

		} else {

			var vertexA = vertices[ vertexIndexA ];
			var vertexB = vertices[ vertexIndexB ];

			edge = {

				a: vertexA, // pointer reference
				b: vertexB,
				newEdge: null,
				// aIndex: a, // numbered reference
				// bIndex: b,
				faces: [] // pointers to face

			};

			map[ key ] = edge;

		}

		edge.faces.push( face );

		metaVertices[ a ].edges.push( edge );
		metaVertices[ b ].edges.push( edge );


	}

	/**
	 * @param {Vector3[]} vertices
	 * @param {Face3[]} faces
	 */
	function generateLookups( vertices, faces, metaVertices, edges ) {

		var i, il, face;

		for ( i = 0, il = vertices.length; i < il; i ++ ) {

			metaVertices[ i ] = { edges: [] };

		}

		for ( i = 0, il = faces.length; i < il; i ++ ) {

			face = faces[ i ];

			processEdge( face.a, face.b, vertices, edges, face, metaVertices );
			processEdge( face.b, face.c, vertices, edges, face, metaVertices );
			processEdge( face.c, face.a, vertices, edges, face, metaVertices );

		}

	}

	function newFace( newFaces, a, b, c, materialIndex ) {

		newFaces.push( new Face3( a, b, c, undefined, undefined, materialIndex ) );

	}

	function midpoint( a, b ) {

		return ( Math.abs( b - a ) / 2 ) + Math.min( a, b );

	}

	function newUv( newUvs, a, b, c ) {

		newUvs.push( [ a.clone(), b.clone(), c.clone() ] );

	}

	/**
 	 * @param {Geometry} geometry
 	 */
	function Loop_Subdivision (geometry){
		var tmp = new Vector3();

		var oldVertices, oldFaces, oldUvs;
		var newVertices, newFaces, newUVs = [];

		var n, i, il, j, k;
		var metaVertices, sourceEdges;

		// new stuff.
		var newEdgeVertices, newSourceVertices;
		
		oldVertices = geometry.vertices; // { x, y, z}
		oldFaces = geometry.faces; // { a: oldVertex1, b: oldVertex2, c: oldVertex3 }
		oldUvs = geometry.faceVertexUvs[ 0 ];

		var hasUvs = oldUvs !== undefined && oldUvs.length > 0;

		/******************************************************
		 *
		 * Step 0: Preprocess Geometry to Generate edges Lookup
		 *
		 *******************************************************/

		metaVertices = new Array( oldVertices.length );
		sourceEdges = {}; // Edge => { oldVertex1, oldVertex2, faces[]  }

		generateLookups( oldVertices, oldFaces, metaVertices, sourceEdges );


		/******************************************************
		 *
		 *	Step 1.
		 *	For each edge, create a new Edge Vertex,
		 *	then position it.
		 *
		 *******************************************************/

		newEdgeVertices = [];
		var other, currentEdge, newEdge, face;
		var edgeVertexWeight, adjacentVertexWeight, connectedFaces;

		for ( i in sourceEdges ) {

			currentEdge = sourceEdges[ i ];
			newEdge = new Vector3();

			edgeVertexWeight = 3 / 8;
			adjacentVertexWeight = 1 / 8;

			connectedFaces = currentEdge.faces.length;

			// check how many linked faces. 2 should be correct.
			if ( connectedFaces != 2 ) {

				// if length is not 2, handle condition
				edgeVertexWeight = 0.5;
				adjacentVertexWeight = 0;

				if ( connectedFaces != 1 ) {

					// console.warn( 'Subdivision Modifier: Number of connected faces != 2, is: ', connectedFaces, currentEdge );

				}

			}

			newEdge.addVectors( currentEdge.a, currentEdge.b ).multiplyScalar( edgeVertexWeight );

			tmp.set( 0, 0, 0 );

			for ( j = 0; j < connectedFaces; j ++ ) {

				face = currentEdge.faces[ j ];

				for ( k = 0; k < 3; k ++ ) {

					other = oldVertices[ face[ ABC[ k ] ] ];
					if ( other !== currentEdge.a && other !== currentEdge.b ) break;

				}

				tmp.add( other );

			}

			tmp.multiplyScalar( adjacentVertexWeight );
			newEdge.add( tmp );

			currentEdge.newEdge = newEdgeVertices.length;
			newEdgeVertices.push( newEdge );

			// console.log(currentEdge, newEdge);

		}

		/******************************************************
		 *
		 *	Step 2.
		 *	Reposition each source vertices.
		 *
		 *******************************************************/

		var beta, sourceVertexWeight, connectingVertexWeight;
		var connectingEdge, connectingEdges, oldVertex, newSourceVertex;
		newSourceVertices = [];

		for ( i = 0, il = oldVertices.length; i < il; i ++ ) {

			oldVertex = oldVertices[ i ];

			// find all connecting edges (using lookupTable)
			connectingEdges = metaVertices[ i ].edges;
			n = connectingEdges.length;

			if ( n == 3 ) {

				beta = 3 / 16;

			} else if ( n > 3 ) {

				beta = 3 / ( 8 * n ); // Warren's modified formula

			}

			// Loop's original beta formula
			// beta = 1 / n * ( 5/8 - Math.pow( 3/8 + 1/4 * Math.cos( 2 * Math. PI / n ), 2) );

			sourceVertexWeight = 1 - n * beta;
			connectingVertexWeight = beta;

			if ( n <= 2 ) {

				// crease and boundary rules
				// console.warn('crease and boundary rules');

				if ( n == 2 ) {

					// console.warn( '2 connecting edges', connectingEdges );
					sourceVertexWeight = 3 / 4;
					connectingVertexWeight = 1 / 8;

					// sourceVertexWeight = 1;
					// connectingVertexWeight = 0;

				} else if ( n == 1 ) {

					// console.warn( 'only 1 connecting edge' );

				} else if ( n == 0 ) {

					// console.warn( '0 connecting edges' );

				}

			}

			newSourceVertex = oldVertex.clone().multiplyScalar( sourceVertexWeight );

			tmp.set( 0, 0, 0 );

			for ( j = 0; j < n; j ++ ) {

				connectingEdge = connectingEdges[ j ];
				other = connectingEdge.a !== oldVertex ? connectingEdge.a : connectingEdge.b;
				tmp.add( other );

			}

			tmp.multiplyScalar( connectingVertexWeight );
			newSourceVertex.add( tmp );

			newSourceVertices.push( newSourceVertex );

		}


		/******************************************************
		 *
		 *	Step 3.
		 *	Generate Faces between source vertices
		 *	and edge vertices.
		 *
		 *******************************************************/

		newVertices = newSourceVertices.concat( newEdgeVertices );
		var sl = newSourceVertices.length, edge1, edge2, edge3;
		newFaces = [];

		var uv, x0, x1, x2;
		var x3 = new Vector2();
		var x4 = new Vector2();
		var x5 = new Vector2();

		for ( i = 0, il = oldFaces.length; i < il; i ++ ) {

			face = oldFaces[ i ];

			// find the 3 new edges vertex of each old face

			edge1 = getEdge( face.a, face.b, sourceEdges ).newEdge + sl;
			edge2 = getEdge( face.b, face.c, sourceEdges ).newEdge + sl;
			edge3 = getEdge( face.c, face.a, sourceEdges ).newEdge + sl;

			// create 4 faces.

			newFace( newFaces, edge1, edge2, edge3, face.materialIndex );
			newFace( newFaces, face.a, edge1, edge3, face.materialIndex );
			newFace( newFaces, face.b, edge2, edge1, face.materialIndex );
			newFace( newFaces, face.c, edge3, edge2, face.materialIndex );

			// create 4 new uv's

			if ( hasUvs ) {

				uv = oldUvs[ i ];

				x0 = uv[ 0 ];
				x1 = uv[ 1 ];
				x2 = uv[ 2 ];

				x3.set( midpoint( x0.x, x1.x ), midpoint( x0.y, x1.y ) );
				x4.set( midpoint( x1.x, x2.x ), midpoint( x1.y, x2.y ) );
				x5.set( midpoint( x0.x, x2.x ), midpoint( x0.y, x2.y ) );

				newUv( newUVs, x3, x4, x5 );
				newUv( newUVs, x0, x3, x5 );

				newUv( newUVs, x1, x4, x3 );
				newUv( newUVs, x2, x5, x4 );

			}

		}

		// Overwrite old arrays
		geometry.vertices = newVertices;
		geometry.faces = newFaces;
		if ( hasUvs ) geometry.faceVertexUvs[ 0 ] = newUVs;

		// console.log('done');
	}

} )();

export { SubdivisionModifier };
