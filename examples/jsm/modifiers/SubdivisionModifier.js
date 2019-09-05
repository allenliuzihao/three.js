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
import { Geometry } from '../../../src/core/Geometry.js';
import { Vector2 } from '../../../src/math/Vector2.js';
import { Vector3 } from '../../../src/math/Vector3.js';


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

	geometry.computeFaceNormals();
	geometry.computeVertexNormals();

	return geometry;

};

( function () {

	/////////////////////////////
	// Performs one iteration of Subdivision
	SubdivisionModifier.prototype.smooth = function ( geometry ) {
		if(this.scheme === "loop") {
			loop (geometry);
		} else if (this.scheme == "catmull"){
			catmull (geometry);
		}
	};

	/* Catmull-Clark subdivision scheme */

	/**
	 * generate "vertex1-vertex2" to faces set per edge
	 * 
	 * @param {Number} vertexAIdx 
	 * @param {Number} vertexBIdx 
	 * @param {Face3} face 
	 * @param {Map<String, Set<Face3>>} edges 
	 */
	function generateEachEdge (vertexAIdx, vertexBIdx, face, edges){
		let vertexIndexA = Math.min( vertexAIdx, vertexBIdx );
		let vertexIndexB = Math.max( vertexAIdx, vertexBIdx );
		let key = vertexIndexA + "_" + vertexIndexB;

		if(!(key in edges)){
			edges.set(key, new Set(face));
		} else {
			edges.get(key).add(face);
		}
	}

	/**
	 * generate "vertex1-vertex2" to faces set per face
	 * 
	 * @param {Face3[]} faces 
	 * @param {Map<String, Set<Face3>>} edges 
	 */
	function generateEdges(faces, edges) {
		let face;

		for (let i = 0, il = faces.length; i < il; i++ ) {
			face = faces[i];
			generateEachEdge(face.a, face.b, face, edges);
			generateEachEdge(face.b, face.c, face, edges);
			generateEachEdge(face.c, face.a, face, edges);
		}
	}

	/**
	 * calculate face point per face
	 * 
	 * @param {Face3[]} faces 
	 * @param {Vector3[]} vertices 
	 * @param {Map<Face3, Vector3>} faceToFacePoints 
	 */
	function generateFacePoints (faces, vertices, faceToFacePoints){
		let face, vertexA, vertexB, vertexC;

		for (let index = 0, length = faces.length; index < length; index++) {
			face = faces[index];
			vertexA = vertices[face.a];
			vertexB = vertices[face.b];
			vertexC = vertices[face.c];
			faceToFacePoints.set(face, 
				new Vector3()
					.addVectors(vertexA, vertexB)
					.add(vertexC)
					.multiplyScalar(1/3));
		}
	}


	function generateEdgePoints (faces , edgeToEdgePoints) {
		
	}

	/**
 	 * @param {Geometry} geometry
 	 */
	function catmull (geometry){
		let faceToFacePoints = new Map(), edgeToEdgePoints = new Map();
		let vertices = geometry.vertices, faces = geometry.faces;
		let edges = new Map();

		// generate edges information; each edge information is "vertex1-vertex2" -> set {faces3}
		generateEdges(faces, edges);

		// generate vertex -> {faces: [], edges: []}

		// calculate the face points for all faces
		generateFacePoints(faces, vertices, faceToFacePoints);

		// calculate edge points
		

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
	function loop (geometry){
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
