using System.Collections.Generic;
using UnityEngine;

namespace InverseKinematicsMesh
{
    public class IK : MonoBehaviour
    {
        public Camera camera;
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
            intersect0.mesh = CreateIntersectRangeMesh(Vector3.zero, shadow0.mesh, shadow0.transform.position, ring0.mesh, ring0.transform.position);
        }

        private void Update()
        {
            //Ray ray = camera.ScreenPointToRay(Input.mousePosition);
            //rayPlane.Raycast(ray, out float hitDist);
            //ring0.transform.position = ray.origin + ray.direction * hitDist;
            if (PointInMesh(target - shadow0.transform.position, shadow0.mesh))
            {
            }
        }

        public bool PointInMesh(Vector3 point, Mesh mesh)
        {
            for (int i = 0; i < mesh.triangles.Length; i += 3)
            {
                Vector3 v0 = mesh.vertices[mesh.triangles[i]];
                Vector3 v1 = mesh.vertices[mesh.triangles[i + 1]];
                Vector3 v2 = mesh.vertices[mesh.triangles[i + 2]];
                if (PointInTriangle(new Vector2(point.z, point.y), 
                    new Vector2(v0.z, v0.y), new Vector2(v1.z, v1.y), new Vector2(v2.z, v2.y)))
                {
                    return true;
                }
            }
            return false;
        }
        public static float sign(Vector2 p1, Vector2 p2, Vector2 p3)
        {
            return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
        }
        public static bool PointInTriangle(Vector2 pt, Vector2 v1, Vector2 v2, Vector2 v3)
        {
            float d1, d2, d3;
            bool has_neg, has_pos;

            d1 = sign(pt, v1, v2);
            d2 = sign(pt, v2, v3);
            d3 = sign(pt, v3, v1);

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

        public Mesh CreateShadowRangeMesh()
        {
            Mesh mesh = new Mesh();

            Vector3[] vertices = new Vector3[GetVertexCount(chain.Length)];
            SetVertex(0, Vector3.zero, vertices);
            // middle line
            float length = 0;
            for (int i = 0; i < chain.Length; i++)
            {
                length += chain[i].length;
                SetVertex(GetVertexIndex(0, i + 1, 0), transform.forward * length, vertices);
            }
            // side lines
            for (int i = 1; i < chain.Length + 1; i++)
            {
                AddLineVertices(i, 0, vertices);
                AddLineVertices(i, 1, vertices);
            }
            mesh.SetVertices(vertices);

            int[] indices = new int[GetTriangleCount(chain.Length) * 3];

            int ti = 0;
            for (int i = 0; i < chain.Length; i++)
            {
                // expanding triangles side 1
                indices[3 * ti]     = GetVertexIndex(i, i, 0);
                indices[3 * ti + 1] = GetVertexIndex(i + 1, i + 1, 0);
                indices[3 * ti + 2] = GetVertexIndex(i, i + 1, 0);
                ti++;

                for (int j = 0; j < i; j++)
                {
                    // center side 1
                    indices[3 * ti] = GetVertexIndex(i - j, i, 0);
                    indices[3 * ti + 1] = GetVertexIndex(i - j, i + 1, 0);
                    indices[3 * ti + 2] = GetVertexIndex(i - j - 1, i + 1, 0);
                    ti++;
                    indices[3 * ti] = GetVertexIndex(i - j, i, 0);
                    indices[3 * ti + 1] = GetVertexIndex(i - j - 1, i + 1, 0);
                    indices[3 * ti + 2] = GetVertexIndex(i - j - 1, i, 0);
                    ti++;
                    // center side 2
                    indices[3 * ti] = GetVertexIndex(i - j, i, 1);
                    indices[3 * ti + 1] = GetVertexIndex(i - j - 1, i + 1, 1);
                    indices[3 * ti + 2] = GetVertexIndex(i - j, i + 1, 1);
                    ti++;
                    indices[3 * ti] = GetVertexIndex(i - j, i, 1);
                    indices[3 * ti + 1] = GetVertexIndex(i - j - 1, i, 1);
                    indices[3 * ti + 2] = GetVertexIndex(i - j - 1, i + 1, 1);
                    ti++;
                }

                // expanding triangles side 2
                indices[3 * ti]     = GetVertexIndex(i, i, 1);
                indices[3 * ti + 1] = GetVertexIndex(i, i + 1, 1);
                indices[3 * ti + 2] = GetVertexIndex(i + 1, i + 1, 1);
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
            vertices[ringVertexCount] = point;
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
        public Mesh CreateIntersectRangeMesh(Vector3 newPos, Mesh shadow, Vector3 shadowPos, Mesh ring, Vector3 ringPos)
        {
            Vector3[] v1 = shadow.vertices; // in world position
            for (int i = 0; i < v1.Length; i++)
            {
                v1[i] += shadowPos;
            }
            float inner = Vector3.Distance(ring.vertices[0], Vector3.zero);
            float outer = Vector3.Distance(ring.vertices[ring.vertexCount - 1], Vector3.zero);
            Vector3[] v2 = ring.vertices; // in world pos
            for (int i = 0; i < v2.Length; i++)
            {
                v2[i] += ringPos;
            }

            Mesh mesh = new Mesh();

            List<Vector3> vertices = new List<Vector3>();
            // vertices of shadow mesh inside ring mesh
            for (int i = 0; i < v1.Length; i++)
            {
                float dist = Vector3.Distance(v1[i], ringPos);
                if (inner <= dist && dist <= outer)
                {
                    vertices.Add(v1[i] - newPos);
                    new GameObject(i.ToString()).transform.position = v1[i] - newPos;
                }
            }
            // vertices of ring mesh inside shadow mesh
            for (int i = 0; i < v2.Length; i++)
            {
                if (PointInMesh(v2[i], shadow))
                {
                    vertices.Add(v2[i] - newPos);
                    new GameObject(i.ToString()).transform.position = v2[i] - newPos;
                }
            }
            // edge outside edge intersection between shadow and ring mesh
            // shadow edges
            List<int> shadowEdges = new List<int>() { 0, 2, 0, 3 };
            for (int i = 1; i < chain.Length; i++)
            {
                shadowEdges.Add(GetVertexIndex(i, i, 0));
                shadowEdges.Add(GetVertexIndex(i + 1, i + 1, 0));

                shadowEdges.Add(GetVertexIndex(i, i, 1));
                shadowEdges.Add(GetVertexIndex(i + 1, i + 1, 1));
            }
            for (int i = 0; i < chain.Length; i++)
            {
                shadowEdges.Add(GetVertexIndex(i, chain.Length, 0));
                shadowEdges.Add(GetVertexIndex(i + 1, chain.Length, 0));

                shadowEdges.Add(GetVertexIndex(i, chain.Length, 1));
                shadowEdges.Add(GetVertexIndex(i + 1, chain.Length, 1));
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
            }


            mesh.SetVertices(vertices.ToArray());

            return mesh;
        }
        private static int GetVertexCount(int arcCount)
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
                return GetVertexCount(arcCount - 1) + 2 * arcCount + 1;
            }
        }
        private static int GetVertexIndex(int startArc, int currentArc, int constraintIndex)
        {
            if (startArc < 0)
            {
                Debug.LogError("Getting vertex index works only for non-negative start arcs");
            }
            if (startArc == 0)
            {
                return GetVertexCount(currentArc - 1);
            }
            switch (constraintIndex)
            {
                case 0:
                    return GetVertexCount(currentArc - 1) + 2 * startArc - 1;
                case 1:
                    return GetVertexCount(currentArc - 1) + 2 * startArc;
                default:
                    Debug.LogError($"Constraint index is {constraintIndex} but only 0 and 1 are valid values");
                    return -1;
            }
        }
        private void AddLineVertices(int startArc, int constraintIndex, Vector3[] vertices)
        {
            Vector3 dir = transform.forward;
            Vector3 startPos = Vector3.zero;
            for (int j = 0; j < startArc; j++)
            {
                dir = Quaternion.AngleAxis(chain[j].constraint[constraintIndex], transform.right) * dir;
                startPos += dir * chain[j].length;
            }

            SetVertex(GetVertexIndex(startArc, startArc, constraintIndex), startPos, vertices);
            float dist = 0;
            for (int i = startArc; i < chain.Length; i++)
            {
                dist += chain[i].length;
                SetVertex(GetVertexIndex(startArc, i + 1, constraintIndex), startPos + dir * dist, vertices);
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