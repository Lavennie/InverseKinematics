using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace InverseKinematicsMesh
{
    public class IK : MonoBehaviour
    {
        const float DEBUG_TIME = 2;

        new public Camera camera;
        public Material shadowMat;
        public Material ringMat;
        public Material intersectMat;

        [Header("IK")]
        public IKBone[] chain;
        public Vector3 target;

        private MeshFilter shadow0;
        private MeshFilter ring0;
        private MeshFilter intersect0;
        private Plane rayPlane;

        private void Awake()
        {
            rayPlane = new Plane(transform.right, Vector3.zero);

            shadow0 = new GameObject("shadow [0]", typeof(MeshFilter), typeof(MeshRenderer)).GetComponent<MeshFilter>();
            shadow0.transform.position = transform.position;
            shadow0.GetComponent<MeshRenderer>().material = shadowMat;
            shadow0.mesh = CreateShadowRangeMesh();

            ring0 = new GameObject("ring [0]", typeof(MeshFilter), typeof(MeshRenderer)).GetComponent<MeshFilter>();
            ring0.GetComponent<MeshRenderer>().material = ringMat;
            ring0.mesh = Create2BoneRangeMesh(chain.Length - 1, 1);
            ring0.transform.position = target;

            intersect0 = new GameObject("intersect [0]", typeof(MeshFilter), typeof(MeshRenderer)).GetComponent<MeshFilter>();
            intersect0.GetComponent<MeshRenderer>().material = intersectMat;
            intersect0.mesh = CreateIntersectRangeMesh(Vector3.zero, shadow0, ring0);
        }

        private void Update()
        {
            //Ray ray = camera.ScreenPointToRay(Input.mousePosition);
            //rayPlane.Raycast(ray, out float hitDist);
            //ring0.transform.position = ray.origin + ray.direction * hitDist;
            if (PointInMesh(target - shadow0.transform.position, shadow0))
            {
            }
        }

        public bool PointInMesh(Vector3 point, MeshFilter mesh, bool debug = false)
        {
            for (int i = 0; i < mesh.mesh.triangles.Length; i += 3)
            {
                Vector3 v0 = mesh.transform.TransformPoint(mesh.mesh.vertices[mesh.mesh.triangles[i]]);
                Vector3 v1 = mesh.transform.TransformPoint(mesh.mesh.vertices[mesh.mesh.triangles[i + 1]]);
                Vector3 v2 = mesh.transform.TransformPoint(mesh.mesh.vertices[mesh.mesh.triangles[i + 2]]);
                if (debug)
                {
                    Debug.DrawLine(v0, v1, Color.green, 1);
                    Debug.DrawLine(v1, v2, Color.green, 1);
                    Debug.DrawLine(v0, v2, Color.green, 1);
                }
                if (PointInTriangle(new Vector2(point.z, point.y), 
                    new Vector2(v0.z, v0.y), new Vector2(v1.z, v1.y), new Vector2(v2.z, v2.y)))
                {
                    return true;
                }
            }
            return false;
        }
        public static float Sign(Vector2 p1, Vector2 p2, Vector2 p3)
        {
            return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
        }
        public static bool PointInTriangle(Vector2 pt, Vector2 v1, Vector2 v2, Vector2 v3)
        {
            float d1, d2, d3;
            bool has_neg, has_pos;

            d1 = Sign(pt, v1, v2);
            d2 = Sign(pt, v2, v3);
            d3 = Sign(pt, v3, v1);

            has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
            has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

            return !(has_neg && has_pos);
        }
        public static bool LineIntersection(Vector2 p1, Vector2 p2, Vector2 p3, Vector2 p4, ref Vector2 intersection)
        {
            float Ax, Bx, Cx, Ay, By, Cy, d, e, f, num/*,offset*/;
            float x1lo, x1hi, y1lo, y1hi;

            Ax = p2.x - p1.x;
            Bx = p3.x - p4.x;

            // X bound box test/
            if (Ax < 0)
            {
                x1lo = p2.x; x1hi = p1.x;
            }
            else
            {
                x1hi = p2.x; x1lo = p1.x;
            }

            if (Bx > 0)
            {
                if (x1hi < p4.x || p3.x < x1lo) return false;
            }
            else
            {
                if (x1hi < p3.x || p4.x < x1lo) return false;
            }

            Ay = p2.y - p1.y;
            By = p3.y - p4.y;

            // Y bound box test//
            if (Ay < 0)
            {
                y1lo = p2.y; y1hi = p1.y;
            }
            else
            {
                y1hi = p2.y; y1lo = p1.y;
            }

            if (By > 0)
            {
                if (y1hi < p4.y || p3.y < y1lo) return false;
            }
            else
            {
                if (y1hi < p3.y || p4.y < y1lo) return false;
            }

            Cx = p1.x - p3.x;
            Cy = p1.y - p3.y;
            d = By * Cx - Bx * Cy;  // alpha numerator//
            f = Ay * Bx - Ax * By;  // both denominator//

            // alpha tests//
            if (f > 0)
            {
                if (d < 0 || d > f) return false;
            }
            else
            {
                if (d > 0 || d < f) return false;
            }

            e = Ax * Cy - Ay * Cx;  // beta numerator//

            // beta tests //
            if (f > 0)
            {
                if (e < 0 || e > f) return false;
            }
            else
            {
                if (e > 0 || e < f) return false;
            }
            // check if they are parallel
            if (f == 0) return false;
            // compute intersection coordinates //
            num = d * Ax; // numerator //
            //    offset = same_sign(num,f) ? f*0.5f : -f*0.5f;   // round direction //
            //    intersection.x = p1.x + (num+offset) / f;
            intersection.x = p1.x + num / f;

            num = d * Ay;
            //    offset = same_sign(num,f) ? f*0.5f : -f*0.5f;
            //    intersection.y = p1.y + (num+offset) / f;
            intersection.y = p1.y + num / f;

            return true;

        }
        public static bool LineMeshClosestIntersection(Vector3 a, Vector3 b, Vector3[] v, int[] e, List<int> ignoreEdges, out Vector3 intersection, out Edge edgeWithIntersect)
        {
            intersection = Vector3.zero;
            edgeWithIntersect = new Edge(-1, Vector3.zero, Vector3.zero);
            List<System.Tuple<Vector3, int>> intersections = new List<System.Tuple<Vector3, int>>();
            for (int j = 0; j < e.Length; j += 2)
            {
                bool ignore = false;
                for (int k = 0; k < ignoreEdges.Count; k += 2)
                {
                    if (e[j] == ignoreEdges[k] && e[j + 1] == ignoreEdges[k + 1])
                    {
                        ignore = true;
                        break;
                    }
                }
                if (ignore) { continue; }

                Vector2 p1 = new Vector2(a.z, a.y);
                Vector2 p2 = new Vector2(b.z, b.y);
                Vector2 p3 = new Vector2(v[e[j]].z, v[e[j]].y);
                Vector2 p4 = new Vector2(v[e[j + 1]].z, v[e[j + 1]].y);
                Vector2 intersect = Vector2.zero;
                if (LineIntersection(p1, p2, p3, p4, ref intersect))
                {
                    intersections.Add(new System.Tuple<Vector3, int>(new Vector3(0, intersect.y, intersect.x), j / 2));
                }
            }
            float minDist = float.MaxValue;
            for (int i = 0; i < intersections.Count; i++)
            {
                float dist = Vector3.Distance(a, intersections[i].Item1);
                if (dist < minDist)
                {
                    minDist = dist;
                    intersection = intersections[i].Item1;
                    edgeWithIntersect = new Edge(intersections[i].Item2, v[e[2 * intersections[i].Item2]], v[e[2 * intersections[i].Item2 + 1]]);
                }
            }
            return intersections.Count > 0;
        }
        public static List<Vector3> LineMeshIntersection(Vector2 p1, Vector2 p2, Vector3[] v, int[] e)
        {
            List<Vector3> intersections = new List<Vector3>();
            for (int j = 0; j < e.Length; j += 2)
            {
                Vector2 p3 = new Vector2(v[e[j]].z, v[e[j]].y);
                Vector2 p4 = new Vector2(v[e[j + 1]].z, v[e[j + 1]].y);
                Vector2 intersection = Vector2.zero;
                if (LineIntersection(p1, p2, p3, p4, ref intersection))
                {
                    intersections.Add(intersection);
                }
            }
            return intersections;
        }

        public Vector3[] GetShadowOutlineVertices(int arcCount, MeshFilter shadow)
        {
            Vector3[] v = new Vector3[GetShadowVertexOutlineCount(arcCount)];
            // origin
            v[0] = shadow.mesh.vertices[0];
            int vi = 1;
            // left side
            for (int i = 1; i < arcCount; i++)
            {
                v[vi] = shadow.mesh.vertices[GetShadowVertexIndex(i, i, 0)];
                vi++;
            }
            // left end
            for (int i = arcCount; i > 0; i--)
            {
                v[vi] = shadow.mesh.vertices[GetShadowVertexIndex(i, arcCount, 0)];
                vi++;
            }
            // center end
            v[vi] = shadow.mesh.vertices[GetShadowVertexIndex(0, arcCount, 0)];
            vi++;
            // right end
            for (int i = 1; i < arcCount + 1; i++)
            {
                v[vi] = shadow.mesh.vertices[GetShadowVertexIndex(i, arcCount, 1)];
                vi++;
            }
            // right side
            for (int i = arcCount - 1; i > 0; i--)
            {
                v[vi] = shadow.mesh.vertices[GetShadowVertexIndex(i, i, 1)];
                vi++;
            }
            return v;
            /*// sides
            for (int i = 0; i < arcCount; i++)
            {
                v[2 * i + 1] = shadow.mesh.vertices[GetShadowVertexIndex(1, 1, 0)];
                v[2 * i + 2] = shadow.mesh.vertices[GetShadowVertexIndex(1, 1, 1)];
            }
            // end
            v[2 * arcCount + 1] = shadow.mesh.vertices[GetShadowVertexIndex(0, arcCount, 0)];
            for (int i = 0; i < arcCount - 1; i++)
            {
                v[2 * arcCount + 2 + 2 * i] = shadow.mesh.vertices[GetShadowVertexIndex(i + 1, arcCount, 0)];
                v[2 * arcCount + 2 + 2 * i + 1] = shadow.mesh.vertices[GetShadowVertexIndex(i + 1, arcCount, 1)];
            }*/
        }

        public void DebugDrawOutline(Vector3[] v, int[] e, Color color)
        {
            for (int i = 0; i < e.Length; i += 2)
            {
                Debug.DrawLine(v[e[i]], v[e[i + 1]], color, DEBUG_TIME);
            }
        }
        public void DebugDrawEdge(Vector3 a, Vector3 b, Color color)
        {
            Debug.DrawLine(a, b, color, DEBUG_TIME);
        }
        public void DebugDrawEdge(Edge edge, Color color)
        {
            Debug.DrawLine(edge.a, edge.b, color, DEBUG_TIME);
        }
        public void DebugDrawEdge(Edge edge, Vector3 origin, Color color)
        {
            Debug.DrawLine(edge.a + origin, edge.b + origin, color, DEBUG_TIME);
        }

        public Mesh CreateShadowRangeMesh()
        {
            Mesh mesh = new Mesh();

            Vector3[] vertices = new Vector3[GetShadowVertexCount(chain.Length)];
            SetVertex(0, Vector3.zero, vertices);
            // middle line
            float length = 0;
            for (int i = 0; i < chain.Length; i++)
            {
                length += chain[i].length;
                SetVertex(GetShadowVertexIndex(0, i + 1, 0), transform.forward * length, vertices);
            }
            // side lines
            for (int i = 1; i < chain.Length + 1; i++)
            {
                AddShadowLineVertices(i, 0, vertices);
                AddShadowLineVertices(i, 1, vertices);
            }
            mesh.SetVertices(vertices);

            int[] indices = new int[GetTriangleCount(chain.Length) * 3];

            int ti = 0;
            for (int i = 0; i < chain.Length; i++)
            {
                // expanding triangles side 1
                indices[3 * ti]     = GetShadowVertexIndex(i, i, 0);
                indices[3 * ti + 1] = GetShadowVertexIndex(i + 1, i + 1, 0);
                indices[3 * ti + 2] = GetShadowVertexIndex(i, i + 1, 0);
                ti++;

                for (int j = 0; j < i; j++)
                {
                    // center side 1
                    indices[3 * ti] = GetShadowVertexIndex(i - j, i, 0);
                    indices[3 * ti + 1] = GetShadowVertexIndex(i - j, i + 1, 0);
                    indices[3 * ti + 2] = GetShadowVertexIndex(i - j - 1, i + 1, 0);
                    ti++;
                    indices[3 * ti] = GetShadowVertexIndex(i - j, i, 0);
                    indices[3 * ti + 1] = GetShadowVertexIndex(i - j - 1, i + 1, 0);
                    indices[3 * ti + 2] = GetShadowVertexIndex(i - j - 1, i, 0);
                    ti++;
                    // center side 2
                    indices[3 * ti] = GetShadowVertexIndex(i - j, i, 1);
                    indices[3 * ti + 1] = GetShadowVertexIndex(i - j - 1, i + 1, 1);
                    indices[3 * ti + 2] = GetShadowVertexIndex(i - j, i + 1, 1);
                    ti++;
                    indices[3 * ti] = GetShadowVertexIndex(i - j, i, 1);
                    indices[3 * ti + 1] = GetShadowVertexIndex(i - j - 1, i, 1);
                    indices[3 * ti + 2] = GetShadowVertexIndex(i - j - 1, i + 1, 1);
                    ti++;
                }

                // expanding triangles side 2
                indices[3 * ti]     = GetShadowVertexIndex(i, i, 1);
                indices[3 * ti + 1] = GetShadowVertexIndex(i, i + 1, 1);
                indices[3 * ti + 2] = GetShadowVertexIndex(i + 1, i + 1, 1);
                ti++;
            }

            mesh.SetTriangles(indices, 0);

            return mesh;
        }
        public Mesh Create2BoneRangeMesh(int endArc, int detail = 3)
        {
            Mesh mesh = new Mesh();
            int ringVertexCount = Mathf.RoundToInt(Mathf.Pow(2, detail + 2));

            Vector3[] vertices = new Vector3[2 * ringVertexCount];

            Quaternion rot = Quaternion.AngleAxis(360.0f / ringVertexCount, transform.right);
            float inner = chain[endArc].length;
            float outer = inner + chain[endArc - 1].length * Mathf.Cos(Mathf.Max(Mathf.Abs
                (chain[endArc - 1].constraint.x), Mathf.Abs(chain[endArc - 1].constraint.y)) * Mathf.Deg2Rad);

            // inner circle
            Vector3 point = transform.forward * inner;
            SetVertex(0, point, vertices);
            for (int i = 1; i < ringVertexCount; i++)
            {
                point = rot * point;
                SetVertex(i, point, vertices);
            }
            // outer circle
            point = transform.forward * outer;
            SetVertex(ringVertexCount, point, vertices);
            for (int i = 1; i < ringVertexCount; i++)
            {
                point = rot * point;
                SetVertex(ringVertexCount + i, point, vertices);
            }

            mesh.SetVertices(vertices);

            int[] indices = new int[3 * vertices.Length];
            int ti = 0;
            for (int i = 0; i < ringVertexCount; i++)
            {
                indices[3 * ti] = i;
                indices[3 * ti + 1] = i + ringVertexCount;
                indices[3 * ti + 2] = (i + 1) % ringVertexCount + ringVertexCount;
                ti++;
                indices[3 * ti] = i;
                indices[3 * ti + 1] = (i + 1) % ringVertexCount + ringVertexCount;
                indices[3 * ti + 2] = (i + 1) % ringVertexCount;
                ti++;
            }

            mesh.SetTriangles(indices, 0);

            return mesh;
        }
        public void FollowOutline(MeshFilter m1, MeshFilter m2, Vector3[] v1, int[] e1, Vector3[] v2, int[] e2, List<Vector3> vertices, List<int> loopEnds)
        {
            // start with edge, check for intersections with mesh2
            // -   there are intersections (do only for closest intersection): continue from edge on which closest intersection is
            // -   there are no intersections: check end vertex of edge, continue to next edge
            //vertices = new List<Vector3>();
            //loopEnds = new List<int>();
            int currentMesh = 0;
            Edge currentEdge;
            List<int> usedEdges = new List<int>();
            HashSet<int> checkedEdges = new HashSet<int>();

            // check all edges of first mesh
            for (int i = 0; i < e1.Length; i += 2)
            {
                if (checkedEdges.Contains(i / 2)) { continue; }
                currentMesh = 0;
                currentEdge = new Edge(i / 2, v1[e1[i]], v1[e1[i + 1]]);
                usedEdges.Clear();
                bool forceAddV = false;
                while (true)
                {
                    /*if (i > -1)
                    {
                        yield return new WaitUntil(() => Input.GetKeyDown(KeyCode.Space));
                        yield return new WaitForSeconds(0.01f);
                    }*/
                    MeshFilter curMesh = (currentMesh == 0) ? m1 : m2;
                    MeshFilter otherMesh = (currentMesh == 0) ? m2 : m1;
                    Vector3[] curV = (currentMesh == 0) ? v1 : v2;
                    Vector3[] otherV = (currentMesh == 0) ? v2 : v1;
                    int[] curE = (currentMesh == 0) ? e1 : e2;
                    int[] otherE = (currentMesh == 0) ? e2 : e1;

                    if (currentMesh == 0)
                    {
                        checkedEdges.Add(currentEdge.index);
                    }

                    bool intersects = LineMeshClosestIntersection(currentEdge.a, currentEdge.b, otherV, otherE, usedEdges, out Vector3 intersection, out Edge intersectEdge);
                    usedEdges = new List<int>() { curE[2 * currentEdge.index], curE[2 * currentEdge.index + 1] };
                    if (intersects)
                    {
                        if (forceAddV)
                        {
                            if (!ContainsApprox(vertices, currentEdge.a))
                            {
                                vertices.Add(currentEdge.a);
                            }
                            else
                            {
                                // a loop was created with outline
                                if (loopEnds.Count == 0 || loopEnds[loopEnds.Count - 1] != vertices.Count)
                                {
                                    loopEnds.Add(vertices.Count);
                                }
                                break;
                            }
                        }
                        ChooseNextOutlineEdge(ref currentEdge, ref currentMesh, usedEdges, intersectEdge, intersection, curMesh, otherMesh, curV, otherV, curE, otherE);
                        forceAddV = true;
                    }
                    else
                    {
                        // whole edge is in other mesh
                        if (PointInMesh(currentEdge.a, otherMesh) || PointInMesh(currentEdge.b, otherMesh))
                        {
                            if (!ContainsApprox(vertices, currentEdge.a))
                            {
                                vertices.Add(currentEdge.a);
                            }
                        }
                        // when whole edge is inside or outside other mesh
                        Edge[] edges = GetOutlineEdgesWithVertex(currentEdge.b, curV, curE);
                        for (int k = 0; k < edges.Length; k++)
                        {
                            if (edges[k].index != currentEdge.index)
                            {
                                if (edges[k].a != currentEdge.b)
                                {
                                    currentEdge = new Edge(edges[k].index, currentEdge.b, edges[k].a);
                                }
                                else
                                {
                                    currentEdge = new Edge(edges[k].index, currentEdge.b, edges[k].b);
                                }
                                break;
                            }
                        }
                        usedEdges.Clear();
                    }
                }
            }
        }
        public bool ContainsApprox(List<Vector3> vertices, Vector3 v)
        {
            const float ERROR = 0.001f; 
            for (int i = 0; i < vertices.Count; i++)
            {
                if (Vector3.Distance(vertices[i], v) < ERROR)
                {
                    return true;
                }
            }
            return false;
        }
        public void ChooseNextOutlineEdge(ref Edge currentEdge, ref int currentMesh, List<int> usedEdges, Edge intersectEdge, Vector3 intersection, MeshFilter curMesh, MeshFilter otherMesh, Vector3[] curV, Vector3[] otherV, int[] curE, int[] otherE)
        {
            bool pointIn1 = PointInMesh(intersectEdge.a, curMesh);
            bool pointIn2 = PointInMesh(intersectEdge.b, curMesh);
            bool intersects1 = LineMeshClosestIntersection(intersection, intersectEdge.a, curV, curE, usedEdges, out Vector3 intersection1, out Edge intersectEdge1);
            bool intersects2 = LineMeshClosestIntersection(intersection, intersectEdge.b, curV, curE, usedEdges, out Vector3 intersection2, out Edge intersectEdge2);

            if (pointIn1 && !intersects1) // choose side 1
            {
                currentEdge = new Edge(intersectEdge.index, intersection, (PointInMesh(intersectEdge.a, curMesh)) ? intersectEdge.a : intersectEdge.b);
            }
            else if (pointIn2 && !intersects2)  // choose side 2
            {
                currentEdge = new Edge(intersectEdge.index, intersection, (PointInMesh(intersectEdge.a, curMesh)) ? intersectEdge.a : intersectEdge.b);
            }
            else if (intersects1)  // choose side 1
            {
                currentEdge = new Edge(intersectEdge.index, intersection, intersection1);
            }
            else if (intersects2) // choose side 2
            {
                currentEdge = new Edge(intersectEdge.index, intersection, intersection2);
            }
            // swap checked mesh
            currentMesh = (currentMesh + 1) % 2;
        }
        public class Edge
        {
            public int index;
            public Vector3 a, b;

            public Edge(int index, Vector3 a, Vector3 b)
            {
                this.index = index;
                this.a = a;
                this.b = b;
            }
        }
        public Edge[] GetOutlineEdgesWithVertex(Vector3 p, Vector3[] v, int[] e)
        {
            int foundCount = 0;
            Edge[] edges = new Edge[2];
            for (int i = 0; i < e.Length; i += 2)
            {
                if (v[e[i]] == p || v[e[i + 1]] == p)
                {
                    edges[foundCount] = new Edge(i / 2, v[e[i]], v[e[i + 1]]);
                    foundCount++;
                    if (foundCount == 2)
                    {
                        return edges;
                    }
                }
            }
            return edges;
        }
        public bool GetFirstUnused(Vector3[] v, List<int> used, out int index)
        {
            for (int i = 0; i < v.Length; i++)
            {
                if (!used.Contains(i))
                {
                    index = i;
                    return true;
                }
            }
            index = -1;
            return false;
        }
        public bool GetFirstUnused(Vector3[] v, List<int> used, out Vector3 vector)
        {
            for (int i = 0; i < v.Length; i++)
            {
                if (!used.Contains(i))
                {
                    vector = v[i];
                    return true;
                }
            }
            vector = Vector3.zero;
            return false;
        }
        public bool GetFirstUnused(int[] e, List<int> used, out int index)
        {
            for (int i = 0; i < e.Length; i++)
            {
                if (!used.Contains(i))
                {
                    index = i;
                    return true;
                }
            }
            index = -1;
            return false;
        }
        public Mesh CreateIntersectRangeMesh(Vector3 newPos, MeshFilter shadow, MeshFilter ring)
        {
            Vector3[] v1 = GetShadowOutlineVertices(chain.Length, shadow); // in world position
            int[] e1 = new int[2 * v1.Length];
            for (int i = 0; i < v1.Length; i++)
            {
                v1[i] = shadow.transform.TransformPoint(v1[i]);
                e1[2 * i] = i;
                e1[2 * i + 1] = (i + 1) % v1.Length;
            }
            //float inner = Vector3.Distance(ring.mesh.vertices[0], Vector3.zero);
            //float outer = Vector3.Distance(ring.mesh.vertices[ring.mesh.vertexCount - 1], Vector3.zero);
            Vector3[] v2 = ring.mesh.vertices; // in world pos
            for (int i = 0; i < v2.Length; i++)
            {
                v2[i] = ring.transform.TransformPoint(v2[i]);
            }
            int circleVertCount = v2.Length / 2;
            int[] e2 = new int[circleVertCount * 4];
            for (int i = 0; i < circleVertCount; i++)
            {
                e2[2 * i] = i;
                e2[2 * i + 1] = (i + 1) % circleVertCount;
            }
            for (int i = 0; i < circleVertCount; i++)
            {
                e2[2 * i + 2 * circleVertCount] = i + circleVertCount;
                e2[2 * i + 1 + 2 * circleVertCount] = (i + 1) % circleVertCount + circleVertCount;
            }

            Mesh mesh = new Mesh();

            List<Vector3> vertices = new List<Vector3>();
            List<int> loopEnds = new List<int>();
            FollowOutline(shadow, ring, v1, e1, v2, e2, vertices, loopEnds);
            List<int> indices = new List<int>();
            TriangulizeOutline(vertices, loopEnds, indices);
            Debug.Log(string.Join(' ', vertices));
            Debug.Log(string.Join(' ', loopEnds));
            /*for (int i = 0; i < loopEnds.Count; i++)
            {
                int start = (i == 0) ? 0 : loopEnds[i - 1];
                int range = loopEnds[i] - start;
                for (int j = 0; j < range; j++)
                {
                    indices.Add(start + j);
                    indices.Add(start + (j + 1) % range);
                }
            }*/
            Debug.Log(string.Join(' ', indices));
            mesh.SetVertices(vertices);
            mesh.SetTriangles(indices, 0);
            //DebugDrawOutline(vertices.ToArray(), indices.ToArray(), Color.green);

            /*Debug.Log(string.Join(' ', vertices));
            Debug.Log(string.Join(' ', loopEnds));
            for (int i = 0; i < loopEnds.Count; i++)
            {
                int start = (i == 0) ? 0 : loopEnds[i];
                for (int k = start; k < loopEnds[i + 1]; k++)
                {
                    Debug.DrawLine(vertices[k], vertices[k + 1], Color.cyan, 2);
                }
            }

            mesh.SetVertices(vertices.ToArray());*/
            return mesh;





            /*// vertices of shadow mesh inside ring mesh
            for (int i = 0; i < v1.Length; i++)
            {
                if (PointInMesh(v1[i], ring, i == 0))
                {
                    vertices.Add(v1[i] - newPos);
                    new GameObject(i.ToString()).transform.position = v1[i] - newPos;
                }
            }
            // vertices of ring mesh inside shadow mesh
            for (int i = 0; i < v2.Length; i++)
            {
                if (PointInMesh(v2[i], shadow, i == 0))
                {
                    vertices.Add(v2[i] - newPos);
                    new GameObject(i.ToString()).transform.position = v2[i] - newPos;
                }
            }
            // edge outside edge intersection between shadow and ring mesh
            // shadow edges
            List<int> shadowEdges = new List<int>();
            for (int i = 0; i < v1.Length; i++)
            {
                shadowEdges.Add(i);
                shadowEdges.Add((i + 1) % v1.Length);
            }
            // ring edges
            List<int> ringEdges = new List<int>();
            int ringVertexCount = v2.Length / 2;
            for (int i = 0; i < ringVertexCount; i++)
            {
                ringEdges.Add(i);
                ringEdges.Add((i + 1) % ringVertexCount);
            }
            for (int i = ringVertexCount; i < 2 * ringVertexCount; i++)
            {
                ringEdges.Add(i);
                ringEdges.Add((i + 1) % ringVertexCount + ringVertexCount);
            }
            // edge intersection
            for (int i = 0; i < shadowEdges.Count; i += 2)
            {
                for (int j = 0; j < ringEdges.Count; j += 2)
                {
                    Vector2 p1 = new Vector2(v1[shadowEdges[i]].z, v1[shadowEdges[i]].y);
                    Vector2 p2 = new Vector2(v1[shadowEdges[i + 1]].z, v1[shadowEdges[i + 1]].y);
                    Vector2 p3 = new Vector2(v2[ringEdges[j]].z, v2[ringEdges[j]].y);
                    Vector2 p4 = new Vector2(v2[ringEdges[j + 1]].z, v2[ringEdges[j + 1]].y);
                    Vector2 intersection = Vector2.zero;
                    if (LineIntersection(p1, p2, p3, p4, ref intersection))
                    {
                        vertices.Add(new Vector3(0, intersection.y, intersection.x) - newPos);
                        new GameObject(i.ToString()).transform.position = vertices[vertices.Count - 1];
                    }
                }
            }*/
        }

        private void TriangulizeOutline(List<Vector3> vertices, List<int> loopEnds, List<int> indices)
        {
            List<int> outlineI = new List<int>();
            for (int i = 0; i < loopEnds.Count; i++)
            {
                int start = (i == 0) ? 0 : loopEnds[i - 1];
                int range = loopEnds[i] - start;
                for (int j = 0; j < range; j++)
                {
                    outlineI.Add(start + j);
                    outlineI.Add(start + (j + 1) % range);
                }
            }
            int[] e = outlineI.ToArray();

            for (int i = 0; i < loopEnds.Count; i++)
            {
                int start = (i == 0) ? 0 : loopEnds[i - 1];
                int range = loopEnds[i] - start;
                int triangleCount = range - 2;
                if (triangleCount == 1)
                {
                    indices.Add(start);
                    indices.Add(start + 1);
                    indices.Add(start + 2);
                    continue;
                }

                Vector3[] v = vertices.GetRange(start, range).ToArray();
                for (int j = 0; j < range; j++)
                {
                    int i0 = start + j;
                    int i1 = start + (j + 1) % range;
                    int i2 = start + (j + 2) % range;
                    Vector2 p0 = new Vector2(vertices[i0].z, vertices[i0].y);
                    Vector2 p2 = new Vector2(vertices[i2].z, vertices[i2].y);

                    List<Vector3> intersections = LineMeshIntersection(p0, p2, v, e);
                    if (intersections.Count == 0)
                    {
                        // check if edge is in/out of outline
                        int i3 = start + (j + 3) % range;
                        Vector3 p3 = new Vector3(vertices[i3].z, vertices[i3].y);
                        Vector3 pi = (p0 + p2) / 2;
                        // TODO: how to check when edge is inside outline
                        // what even is inside of outline? define when creating the outline
                        if (LineMeshIntersection(pi, p3, v, e).Count == 0)
                        {

                        }
                    }

                    indices.Add(start);
                    indices.Add(start + j + 1);
                    indices.Add(start + j + 2);
                }
            }
        }
        private static int GetShadowVertexOutlineCount(int arcCount)
        {
            return 4 * arcCount;
        }
        private static int GetShadowVertexCount(int arcCount)
        {
            if (arcCount < 0)
            {
                return 0;
            }
            if (arcCount == 0)
            {
                return 1;
            }
            else
            {
                return GetShadowVertexCount(arcCount - 1) + 2 * arcCount + 1;
            }
        }
        private static int GetShadowVertexIndex(int startArc, int currentArc, int constraintIndex)
        {
            if (startArc < 0)
            {
                Debug.LogError("Getting vertex index works only for non-negative start arcs");
            }
            if (startArc == 0)
            {
                return GetShadowVertexCount(currentArc - 1);
            }
            switch (constraintIndex)
            {
                case 0:
                    return GetShadowVertexCount(currentArc - 1) + 2 * startArc - 1;
                case 1:
                    return GetShadowVertexCount(currentArc - 1) + 2 * startArc;
                default:
                    Debug.LogError($"Constraint index is {constraintIndex} but only 0 and 1 are valid values");
                    return -1;
            }
        }
        private void AddShadowLineVertices(int startArc, int constraintIndex, Vector3[] vertices)
        {
            Vector3 dir = transform.forward;
            Vector3 startPos = Vector3.zero;
            for (int j = 0; j < startArc; j++)
            {
                dir = Quaternion.AngleAxis(chain[j].constraint[constraintIndex], transform.right) * dir;
                startPos += dir * chain[j].length;
            }

            SetVertex(GetShadowVertexIndex(startArc, startArc, constraintIndex), startPos, vertices);
            float dist = 0;
            for (int i = startArc; i < chain.Length; i++)
            {
                dist += chain[i].length;
                SetVertex(GetShadowVertexIndex(startArc, i + 1, constraintIndex), startPos + dir * dist, vertices);
            }
        }
        private void SetVertex(int index, Vector3 vec, Vector3[] vertices)
        {
            new GameObject(index.ToString()).transform.position = transform.position + vec;
            vertices[index] = vec;
        }

        private static int GetTriangleCount(int arcCount)
        {
            if (arcCount <= 0)
            {
                return 0;
            }
            else
            {
                return GetTriangleCount(arcCount - 1) + arcCount * 4 - 2;
            }
        }
    }
}