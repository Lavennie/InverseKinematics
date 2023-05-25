using System.Collections.Generic;
using UnityEditor;
using UnityEngine;

namespace InverseKinematics2D
{
    /// <summary>
    /// 2D IK in local space defined by forward (x) and up (y) axises. Rotations are made about the right axis (z)
    /// </summary>
    public class IK2 : MonoBehaviour
    {
        public Transform target;
        public IK2Segment[] segments;

        private Reach[] dataFront;
        private Reach[] dataBack;

        private Vector2 mousePos = Vector2.zero;
        private Vector2 reach = Vector2.zero;

        private void Awake()
        {
            Init();
        }
        private void Init()
        {
            dataFront = new Reach[segments.Length];
            int i = 0;
            for (int k = segments.Length; k > 0; k--)
            {
                dataFront[i] = new Reach(this, new IndexRange(0, k));
                i++;
            }
            dataBack = new Reach[segments.Length];
            for (int j = 0; j < segments.Length; j++)
            {
                dataBack[j] = new Reach(this, new IndexRange(j, segments.Length));
            }
        }

#if UNITY_EDITOR
        private void OnGUI()
        {
            EditorGUI.LabelField(new Rect(0, 0, 300, EditorGUIUtility.singleLineHeight),
                $"Mouse position: {mousePos}");
            EditorGUI.LabelField(new Rect(0, EditorGUIUtility.singleLineHeight, 300, EditorGUIUtility.singleLineHeight),
                $"Reach: {reach}");
        }
#endif
        private void Update()
        {
            if (target == null)
            {
                mousePos = Camera.main.ScreenToWorldPoint(new Vector3(Input.mousePosition.x, Input.mousePosition.y, -Camera.main.transform.position.z));
            }
            else
            {
                mousePos = target.localPosition;
            }
            mousePos = ClampTarget(mousePos);
            float angle = Mathf.Asin(mousePos.y / mousePos.magnitude) * Mathf.Rad2Deg;

            Interval reachInterval = dataFront[0].Get(this, angle);
            reach = new Vector2(reachInterval.min, reachInterval.max);
            Debug.DrawRay(transform.position + transform.TransformVector(mousePos.normalized * reach.x), transform.TransformVector(mousePos.normalized * Mathf.Max(reach.y - reach.x, Mathf.Epsilon)), Color.red, 0.01f);

            Align(Vector2.zero, mousePos, -45, segments.Length);
        }

        public void Align(Vector3 origin, Vector2 target, float angle, int n, float prevAngle = 0)
        {
            if (n == 1)
            {
                Debug.DrawRay(transform.TransformPoint(origin), transform.TransformVector((target - (Vector2)origin).normalized * segments[segments.Length - 1].length), Color.cyan);
                return;
            }
            Debug.DrawPoint(transform.TransformPoint(target), Color.white, 3);

            int segmentI = segments.Length - n;
            Circle targetCircle = new Circle(new Vector2(-segments[segmentI].length, 0), Vector3.Distance(origin, target));

            Reach smallerReach = dataBack[segmentI + 1];
            float amin = segments[segmentI].angleLimit.min;
            float amax = segments[segmentI].angleLimit.max;
            Vector3 c1 = Vector(amin, segments[segmentI].length);
            Vector3 c2 = Vector(amax, segments[segmentI].length);
            float rotDebugAngle = prevAngle + ((segments[segmentI].targetAngle * 2) - 1) * 180;
            Vector2 localTarget = new Vector2(float.NaN, float.NaN);
            Vector2 toTarget = target - (Vector2)origin;
            float ta = Mathf.DeltaAngle(prevAngle, Angle(toTarget));
            localTarget = Vector(ta, toTarget.magnitude);

            // get intersection intervals between target circle and reach (facing directly right)
            List<Vector2> minPoints = new List<Vector2>();
            List<Vector2> maxPoints = new List<Vector2>();
            foreach (var circleInterval in smallerReach.GetMinCircles(this))
            {
                Circle circle = circleInterval.Key;
                CircleIntersection Inter = circle.IntersectWithCircle(targetCircle);

                switch (Inter.Type)
                {
                    case CircleIntersection.Variant.Touching:
                        minPoints.Add(Inter.p1 + new Vector2(segments[segmentI].length, 0));
                        break;
                    case CircleIntersection.Variant.Intersect:
                        float a1 = Angle(Inter.p1 - circle.center);
                        float a2 = Angle(Inter.p2 - circle.center);
                        if (circleInterval.Value.Contains(a1))
                        {
                            minPoints.Add(Inter.p1);
                        }
                        if (circleInterval.Value.Contains(a2))
                        {
                            minPoints.Add(Inter.p2);
                        }
                        break;
                }
            }
            foreach (var circleInterval in smallerReach.GetMaxCircles(this))
            {
                Circle circle = circleInterval.Key;
                CircleIntersection Inter = circle.IntersectWithCircle(targetCircle);

                switch (Inter.Type)
                {
                    case CircleIntersection.Variant.Touching:
                        maxPoints.Add(Inter.p1 + new Vector2(segments[segmentI].length, 0));
                        break;
                    case CircleIntersection.Variant.Intersect:
                        float a1 = Angle(Inter.p1 - circle.center);
                        float a2 = Angle(Inter.p2 - circle.center);
                        if (circleInterval.Value.Contains(a1))
                        {
                            maxPoints.Add(Inter.p1);
                        }
                        if (circleInterval.Value.Contains(a2))
                        {
                            maxPoints.Add(Inter.p2);
                        }
                        break;
                }
            }

            #region multiinterval magic

            if (maxPoints.Count == 0)
            {
                float a = Angle(target);
                a = Interval.ClampTo(a, smallerReach.ValidInterval);
                maxPoints.Add(Vector(a, smallerReach.Get(this, a).max));
            }
            else if (maxPoints.Count == 1)
            {
                Vector2 v = FollowChain(new IndexRange(0, n - 1), IntervalSide.Min, out _);
                if (Vector2.Distance(targetCircle.center, v) <= targetCircle.radius)
                {
                    maxPoints.Add(v);
                }
                v = FollowChain(new IndexRange(0, n - 1), IntervalSide.Max, out _);
                if (Vector2.Distance(targetCircle.center, v) <= targetCircle.radius)
                {
                    maxPoints.Add(v);
                }
            }
            else
            {
                if (Angle(maxPoints[0]) > Angle(maxPoints[1]))
                {
                    Vector2 temp = maxPoints[0];
                    maxPoints[0] = maxPoints[1];
                    maxPoints[1] = temp;
                }
            }
            minPoints.Sort((a, b) => Angle(a).CompareTo(Angle(b)));

            Interval newTargetInterval = new Interval(Angle(maxPoints[0]), Angle(maxPoints[1 % maxPoints.Count]));

            List<Interval> intervals = new List<Interval>();
            if (maxPoints.Count == 1)
            {
                if (minPoints.Count == 0)
                {
                    intervals.Add(new Interval(Angle(maxPoints[0] - targetCircle.center), Angle(maxPoints[0] - targetCircle.center)));
                }
                else
                {
                    intervals.Add(new Interval(Angle(minPoints[0] - targetCircle.center), Angle(maxPoints[0] - targetCircle.center)));
                }
            }
            else
            {
                if (minPoints.Count == 0)
                {
                    intervals.Add(new Interval(Angle(maxPoints[0] - targetCircle.center), Angle(maxPoints[1] - targetCircle.center)));
                }
                else if (minPoints.Count % 2 == 0)
                {
                    intervals.Add(new Interval(Angle(maxPoints[0] - targetCircle.center), Angle(minPoints[0] - targetCircle.center)));
                    for (int i = 1; i < minPoints.Count - 1; i += 2)
                    {
                        intervals.Add(new Interval(Angle(minPoints[i] - targetCircle.center), Angle(minPoints[i + 1] - targetCircle.center)));
                    }
                    intervals.Add(new Interval(Angle(minPoints[minPoints.Count - 1] - targetCircle.center), Angle(maxPoints[1] - targetCircle.center)));
                }
                else
                {
                    if (maxPoints[0].magnitude < maxPoints[1].magnitude)
                    {
                        for (int i = 0; i < minPoints.Count - 1; i += 2)
                        {
                            intervals.Add(new Interval(Angle(minPoints[i]), Angle(minPoints[i + 1])));
                        }
                        intervals.Add(new Interval(Angle(minPoints[minPoints.Count - 1]), Angle(maxPoints[1])));
                    }
                    else
                    {
                        intervals.Add(new Interval(Angle(maxPoints[0]), Angle(minPoints[0])));
                        for (int i = 0; i < minPoints.Count - 1; i += 2)
                        {
                            intervals.Add(new Interval(Angle(minPoints[i]), Angle(minPoints[i + 1])));
                        }
                    }
                }
            }
            #endregion

            float targetAngle;
            if (segmentI > 0)
            {
                targetAngle = Angle(target - (Vector2)(origin - (Vector3)Vector(prevAngle, segments[segmentI - 1].length)));
            }
            else
            {
                targetAngle = Angle(target);
            }
            targetAngle = Angle(localTarget);
            float partAngle;
            if (intervals.Count == 1 && intervals[0].IsPoint())
            {
                partAngle = prevAngle + targetAngle - intervals[0].min;
            }
            else
            {
                List<Interval> angleOffsets = new List<Interval>();
                for (int i = 0; i < intervals.Count; i++)
                {
                    Interval valid = Interval.Overlap(new Interval(targetAngle - intervals[i].min, targetAngle - intervals[i].max), segments[segmentI].angleLimit);

                    if (!valid.IsEmpty())
                    {
                        angleOffsets.Add(valid);
                    }
                }
                angleOffsets.Sort((a, b) => a.min.CompareTo(b.min));

                MultiInterval mi = new MultiInterval(angleOffsets);
                partAngle = prevAngle + mi.GetConnected(segments[segmentI].targetAngle);
            }
            Vector3 partVec = Vector(partAngle, segments[segmentI].length);
            Debug.DrawRay(transform.TransformPoint(origin), transform.TransformVector(partVec), Color.cyan);

            Align(origin + partVec, target, angle, n - 1, partAngle);
        }

        private Vector2 ClampTarget(Vector2 target)
        {
            float angle = Angle(target);
            float radius = target.magnitude;

            angle = Interval.ClampTo(angle, dataFront[0].ValidInterval);
            radius = Interval.ClampTo(radius, dataFront[0].Get(this, angle));
            return radius * new Vector2(Mathf.Cos(angle * Mathf.Deg2Rad), Mathf.Sin(angle * Mathf.Deg2Rad));
        }

        /// <summary>
        /// Get the local end point if first couple of segments are completely bend to one side and the rest to the other side
        /// <para>[<see cref="IndexRange.startN"/>, <see cref="IndexRangeK.k"/>) -> <paramref name="limitIndex1"/></para>
        /// <para>[<see cref="IndexRangeK.k"/>, <see cref="IndexRange.endN"/>) -> <paramref name="limitIndex2"/></para>
        /// </summary>
        /// <param name="fromToSplit">Bend indices from <see cref="IndexRange.startN"/> (included) to <see cref="IndexRangeK.k"/> (excluded) to side given with <paramref name="limitIndex1"/> 
        /// and from <see cref="IndexRangeK.k"/> (included) to <see cref="IndexRange.endN"/> (excluded) to side given with <paramref name="limitIndex2"/></param>
        /// <param name="limitIndex1">Side to which to bend segments from <see cref="IndexRange.startN"/> (included) to <see cref="IndexRangeK.k"/> (excluded)</param>
        /// <param name="limitIndex2">Side to which to bend <see cref="IndexRangeK.k"/> (included) and all segments that come after until <see cref="IndexRange.endN"/> (excluded)</param>
        /// <returns>Point reached by chain after bending relative to parent of first bend segment (or root)</returns>
        public Vector2 FollowChain(IndexRangeK fromToSplit, IntervalSide limitIndex1, IntervalSide limitIndex2)
        {
            Vector2 p1 = FollowChain(new IndexRange(fromToSplit.startN, fromToSplit.k), limitIndex1, out float endAngle1);
            Vector2 p2 = FollowChain(new IndexRange(fromToSplit.k, fromToSplit.endN), limitIndex2, out _);
            // rotate around normal z axis because p1 and p2 are local vectors
            return p1 + (Vector2)(Quaternion.Euler(0, 0, endAngle1) * p2);
        }
        /// <summary>
        /// Bend subset of segments from <see cref="IndexRange.startN"/> (included) to <see cref="IndexRange.endN"/> (excluded) to side given with <paramref name="limitIndex"/>
        /// <para>[<see cref="IndexRange.startN"/>, <see cref="IndexRange.endN"/>) -> <paramref name="limitIndex"/></para>
        /// </summary>
        /// <param name="fromTo">Bend indices between <see cref="IndexRange.startN"/> (included) and <see cref="IndexRange.endN"/> (excluded)</param>
        /// <param name="limitIndex">Side to which to bend the segments</param>
        /// <param name="endAngle">Angle of forward vector of the last bend segment relative to parent of first bend segment (or root)</param>
        /// <returns>Point reached by chain after bending relative to parent of first bend segment (or root)</returns>
        public Vector2 FollowChain(IndexRange fromTo, IntervalSide limitIndex, out float endAngle)
        {
            endAngle = 0;
            Vector2 p = Vector2.zero;
            for (int i = fromTo.startN; i < fromTo.endN; i++)
            {
                if (limitIndex == IntervalSide.Min || limitIndex == IntervalSide.Max)
                {
                    endAngle += segments[i].angleLimit[(int)limitIndex];
                }
                p += Vector(endAngle, segments[i].length);
            }
            return p;
        }

        /// <summary>
        /// Create a vector in IK space at angle with length
        /// </summary>
        /// <param name="angle">Angle of vector in IK space</param>
        /// <param name="length">Length of vector</param>
        /// <returns>A vector in IK space</returns>
        public static Vector2 Vector(float angle, float length)
        {
            return new Vector2(Mathf.Cos(angle * Mathf.Deg2Rad), Mathf.Sin(angle * Mathf.Deg2Rad)) * length;
        }
        /// <summary>
        /// Calculates the angle to a point in IK space
        /// </summary>
        /// <param name="p">Point given in the IK space</param>
        /// <returns>Angle to point in IK space</returns>
        public static float Angle(Vector2 p)
        {
            return Vector2.SignedAngle(Vector2.right, p);
        }
        /// <summary>
        /// Calculate the angle interval between two points in IK space
        /// </summary>
        /// <param name="p1">Point in IK space</param>
        /// <param name="p2">Point in IK space</param>
        /// <returns>Angle interval between <paramref name="p1"/> and <paramref name="p2"/></returns>
        public static Interval AngleInterval(Vector2 p1, Vector2 p2)
        {
            return new Interval(Angle(p1), Angle(p2));
        }

        /// <summary>
        /// Number of segments in this IK chain
        /// </summary>
        public int SegmentCount { get { return segments.Length; } }
        /// <summary>
        /// Get segment by index
        /// </summary>
        /// <param name="i">Index to get segment by</param>
        /// <returns>Segment at index</returns>
        public IK2Segment this[int i] { get { return segments[i]; } }
        /// <summary>
        /// IK space x axis vector in 3D
        /// </summary>
        public Vector3 Right { get { return transform.forward; } }
        /// <summary>
        /// IK space y axis vector in 3D
        /// </summary>
        public Vector3 Up { get { return transform.up; } }
        /// <summary>
        /// IK space z axis vector in 3D
        /// </summary>
        public Vector3 RotateAxis { get { return transform.right; } }
    }

    /// <summary>
    /// A range of indices given with start and end index
    /// </summary>
    public class IndexRange
    {
        /// <summary>
        /// Index that starts the range
        /// </summary>
        public int startN;
        /// <summary>
        /// Index that ends the range
        /// </summary>
        public int endN;

        public IndexRange(int startN, int partCount)
        {
            this.startN = startN;
            this.endN = partCount;
        }

        /// <summary>
        /// <see cref="endN"/> - <see cref="startN"/>
        /// </summary>
        public int Size { get { return endN - startN; } }
    }
    /// <summary>
    /// Range of indices made up of two ranges that are connected with a divider index somewhere inbetween
    /// </summary>
    public class IndexRangeK : IndexRange
    {
        /// <summary>
        /// The index that divides the two ranges
        /// </summary>
        public int k;

        public IndexRangeK(int startN, int k, int endN) : base(startN, endN)
        {
            this.k = k;
        }
    }
    public enum IntervalSide : byte { Min = 0, Max = 1, Center = 2 }
    [System.Serializable]
    public struct Interval
    {
        public float min;
        public float max;

        public Interval(float min, float max)
        {
            if (min < max)
            {
                this.min = min;
                this.max = max;
            }
            else
            {
                this.min = max;
                this.max = min;
            }
        }

        public static bool AreOverlapping(Interval i1, Interval i2)
        {
            // one end of i2 is inside i1
            return !Overlap(i1, i2).IsEmpty();
            //return i1.min == i2.min && i1.max == i2.max || i1.min < i2.min && i2.min < i1.max || i1.min < i2.max && i2.max < i1.max;
        }
        public static float ClampTo(float value, Interval interval)
        {
            return Mathf.Clamp(value, interval.min, interval.max);
        }
        public static Interval Overlap(Interval i1, Interval i2)
        {
            if (i1.min == i2.min && i1.max == i2.max)
            {
                return i1;
            }
            else if (i2.max < i1.min || i1.max < i2.min)
            {
                return Interval.Empty;
            }
            else
            {
                return new Interval(Mathf.Max(i1.min, i2.min), Mathf.Min(i1.max, i2.max));
            }
        }

        public bool IsEmpty()
        {
            return float.IsNaN(min) && float.IsNaN(max);
        }
        public bool IsPoint()
        {
            return min == max;
        }
        public bool Contains(float angle)
        {
            return min <= angle && angle <= max;
        }
        public float RandomInside()
        {
            return Random.Range(min, max);
        }

        public override string ToString()
        {
            return $"[{min}, {max}]";
        }

        public static Interval operator +(Interval i, float a)
        {
            return new Interval(i.min + a, i.max + a);
        }
        public static Interval operator +(float a, Interval i)
        {
            return new Interval(i.min + a, i.max + a);
        }

        public static Interval Empty { get { return new Interval(float.NaN, float.NaN); } }
        public float this[int i] { get { return (i % 2 == 0) ? min : max; } }
    }
    public class MultiInterval
    {
        private List<Interval> intervals;
        private float sum;

        public MultiInterval(IEnumerable<Interval> intervals)
        {
            this.intervals = new List<Interval>(intervals);
            sum = 0;
            foreach (var i in intervals)
            {
                sum += i.max - i.min;
            }
        }
        public MultiInterval(params Interval[] intervals) : this((IEnumerable<Interval>)intervals) { }

        public float GetConnected(float t)
        {
            float value = float.NaN;
            float tempValue = Mathf.Lerp(0, 1, t);
            float prevEnd = 0;
            for (int i = 0; i < intervals.Count; i++)
            {
                Interval tempInterval = new Interval(prevEnd, prevEnd + (intervals[i].max - intervals[i].min) / sum);
                if (i == intervals.Count - 1)
                {
                    tempInterval.max = 1.0f;
                }
                prevEnd = tempInterval.max;

                if (tempInterval.Contains(tempValue))
                {
                    tempValue = (tempValue - tempInterval.min) / (tempInterval.max - tempInterval.min);
                    value = Mathf.Lerp(intervals[i].min, intervals[i].max, tempValue);
                    break;
                }
            }
            return value;
        }
    }
    /// <summary>
    /// A 2D circle represented with center and radius
    /// </summary>
    public struct Circle
    {
        public Vector2 center;
        public float radius;

        public Circle(Vector2 center, float radius)
        {
            this.center = center;
            this.radius = radius;
        }

        /// <summary>
        /// Intersection with line that goes through origin and has a specific angle
        /// </summary>
        /// <param name="angle">Angle of line</param>
        /// <returns>Two ntersection points with line saved in <see cref="CircleIntersection"/> that is always <see cref="CircleIntersection.Variant.Intersect"/></returns>
        public CircleIntersection IntersectWithLine(float angle)
        {
            float radAngle = angle * Mathf.Deg2Rad;
            float cos = Mathf.Cos(radAngle);
            float sin = Mathf.Sin(radAngle);
            float y1 = sin * (center.y * sin + center.x * cos + Mathf.Sqrt(Mathf.Pow(center.y * sin + center.x * cos, 2) - center.x * center.x - center.y * center.y + radius * radius));
            float y2 = sin * (center.y * sin + center.x * cos - Mathf.Sqrt(Mathf.Pow(center.y * sin + center.x * cos, 2) - center.x * center.x - center.y * center.y + radius * radius));
            // TODO: sign depending on angle
            float x1 = center.x + Mathf.Sqrt(radius * radius - (y1 - center.y) * (y1 - center.y));
            float x2 = center.x - Mathf.Sqrt(radius * radius - (y2 - center.y) * (y2 - center.y));

            return new CircleIntersection(x1, y1, x2, y2);
        }
        public CircleIntersection IntersectWithCircle(Circle circle)
        {
            Vector2 c1 = center;
            Vector2 c2 = circle.center;
            float r1 = radius;
            float r2 = circle.radius;
            if (Vector2.Distance(c1, c2) > r1 + r2 + 2 * Mathf.Epsilon)
            {
                return CircleIntersection.Miss;
            }
            if (c1.y != c2.y)
            {
                float k = -(c1.x - c2.x) / (c1.y - c2.y);
                float a = -((r1 * r1 - r2 * r2) - (c1.x * c1.x - c2.x * c2.x) - (c1.y * c1.y - c2.y * c2.y)) / (2 * (c1.y - c2.y));

                float d = ((c1.x + k * c1.y - k * a) * (c1.x + k * c1.y - k * a)) -
                    (k * k + 1) * (c1.x * c1.x + c1.y * c1.y + a * a - 2 * a * c1.y - r1 * r1);
                // quadratic formula
                float x1 = ((c1.x + k * c1.y - k * a) + Mathf.Sqrt(d)) / (k * k + 1);
                float x2 = ((c1.x + k * c1.y - k * a) - Mathf.Sqrt(d)) / (k * k + 1);
                float y1 = k * x1 + a;
                float y2 = k * x2 + a;
                return new CircleIntersection(x1, y1, x2, y2);
            }
            else if (c1.x != c2.x)
            {
                // law of cosines without acos, since it gets canceled out when used in next line
                float d = Mathf.Abs(c1.x - c2.x);
                if (Mathf.Abs(r1 + r2 - d) < Mathf.Epsilon)
                {
                    return new CircleIntersection((c1.x < c2.x) ? c1.x + r1 : c1.x - r1, c1.y, float.NaN, float.NaN);
                }
                else if (Mathf.Abs(r1 + d - r2) < Mathf.Epsilon)
                {
                    return new CircleIntersection((c1.x < c2.x) ? c1.x - r1 : c1.x + r1, c1.y, float.NaN, float.NaN);
                }
                else if (Mathf.Abs(r2 + d - r1) < Mathf.Epsilon)
                {
                    return new CircleIntersection((c1.x < c2.x) ? c1.x + r1 : c1.x - r1, c1.y, float.NaN, float.NaN);
                }
                float angle = (r1 * r1 + d * d - r2 * r2) / (2 * r1 * d);
                float offset = r1 * angle;
                float x = (c1.x < c2.x) ? c1.x + offset : c1.x - offset;
                CircleIntersection intersection = new Circle(new Vector2(x - c1.x, 0), r1).IntersectWithLine(90);
                intersection.p1 = new Vector2(x, intersection.p1.y + c1.y);
                intersection.p2 = new Vector2(x, intersection.p2.y + c1.y);
                return intersection;
            }
            else
            {
                return (r1 == r2) ? new CircleIntersection(c1.x, c1.y, r1, float.PositiveInfinity) : CircleIntersection.Miss;
            }
        }

        public override string ToString()
        {
            return $"({center}, {radius})";
        }
    }
    /// <summary>
    /// Data about a circle to circle intersection
    /// </summary>
    public struct CircleIntersection
    {
        /// <summary>
        /// Type of circle intersection by number of intersection points 
        /// </summary>
        public enum Variant { 
            /// <summary>
            /// Circles do not intersect
            /// </summary>
            Miss, 
            /// <summary>
            /// The circles are touching in exactly one point
            /// </summary>
            Touching, 
            /// <summary>
            /// Circles intersect in two points
            /// </summary>
            Intersect, 
            /// <summary>
            /// Circles perfectly cover each other (are the same circle), they intersect in all points on edge
            /// </summary>
            Covering 
        }

        public Vector2 p1, p2;

        public CircleIntersection(float x1, float y1, float x2, float y2) : this(new Vector2(x1, y1), new Vector2(x2, y2)) { }
        public CircleIntersection(Vector2 p1, Vector2 p2)
        {
            this.p1 = p1;
            this.p2 = p2;
        }

        /// <summary>
        /// Type of intersection. Defined by whether <see cref="p1"/> or/and <see cref="p2"/> are <see cref="float.NaN"/> or <see cref="float.PositiveInfinity"/>.
        /// <para>If p1.x == <see cref="float.NaN"/> and p2.x == <see cref="float.NaN"/> then it is <see cref="Variant.Miss"/></para>
        /// <para>If p1.x == <see cref="float.NaN"/> and p2.x != <see cref="float.NaN"/> then it is <see cref="Variant.Touching"/>, point is stored in <see cref="p1"/></para>
        /// <para>If p1.x != <see cref="float.NaN"/> and p2.y == <see cref="float.PositiveInfinity"/> then it is <see cref="Variant.Covering"/>, get values from circle and radius</para>
        /// <para>If p1.x != <see cref="float.NaN"/> and p2.y != <see cref="float.PositiveInfinity"/> then it is <see cref="Variant.Intersect"/>, points are stored in <see cref="p1"/> and <see cref="p2"/></para>
        /// </summary>
        public Variant Type
        {
            get
            {
                if (float.IsNaN(p2.x))
                {
                    if (float.IsNaN(p1.x))
                    {
                        return Variant.Miss;
                    }
                    else
                    {
                        return Variant.Touching;
                    }
                }
                else if (float.IsInfinity(p2.y))
                {
                    return Variant.Covering;
                }
                else
                {
                    return Variant.Intersect;
                }
            }
        }
        /// <summary>
        /// Create instance that represents a missing intersection
        /// </summary>
        public static CircleIntersection Miss { get { return new CircleIntersection(float.NaN, float.NaN, float.NaN, float.NaN); } }
    }



    [System.Serializable]
    public class Reach
    {
        private MinReachData minData;
        private MaxReachData maxData;

        public Reach(IK2 ik, IndexRange segmentRange)
        {
            minData = new MinReachData(ik, segmentRange);
            maxData = new MaxReachData(ik, segmentRange);
        }

        public IEnumerable<KeyValuePair<Circle, Interval>> GetMinCircles(IK2 ik)
        {
            for (int i = 0; i < minData.AngleCount - 1; i++)
            {
                yield return new KeyValuePair<Circle, Interval>(minData.GetCircle(ik, i), minData.GetCircleInterval(ik, i));
            }
        }
        public IEnumerable<KeyValuePair<Circle, Interval>> GetMaxCircles(IK2 ik)
        {
            for (int i = -(maxData.AngleCount - 1); i < maxData.AngleCount; i++)
            {
                yield return new KeyValuePair<Circle, Interval>(maxData.GetCircle(ik, i), maxData.GetCircleInterval(ik, i));
            }
        }

        public Interval Get(IK2 ik, float angle)
        {
            return new Interval(minData.Get(ik, angle), maxData.Get(ik, angle));
        }
        public Interval ValidInterval { get { return maxData.ValidInterval; } }
    }

    /// <summary>
    /// Precalculated data from which min reach at an angle is calculated
    /// </summary>
    [System.Serializable]
    public class MinReachData
    {
        private float[] angles;
        private MinReachRegion[] regions;
        private IndexRange range;

        public MinReachData(IK2 ik, IndexRange range)
        {
            this.range = new IndexRange(range.startN, range.endN);

            if (range.Size <= 0) { return; }
            // only one chain part
            if (range.Size == 1)
            {
                this.angles = new float[2] { ik[range.startN].angleLimit.min, ik[range.startN].angleLimit.max };
                this.regions = new MinReachRegion[1] { new MinReachRegion(IntervalSide.Min, range.startN) };
                return;
            }
            // [0]: shorter bend side
            // [1]: longer bend side
            IntervalSide shorterSide = IntervalSide.Min;
            Vector3[] bendPoints = new Vector3[2]
            {
                ik.FollowChain(range, IntervalSide.Min, out _),
                ik.FollowChain(range, IntervalSide.Max, out _)
            };
            float[] lengths = new float[2]
            {
                bendPoints[0].magnitude,
                bendPoints[1].magnitude
            };
            List<float>[] angles = new List<float>[2] { new List<float>(), new List<float>() };
            List<MinReachRegion>[] regions = new List<MinReachRegion>[] { new List<MinReachRegion>(), new List<MinReachRegion>() };
            int[] insertI = new int[2];

            // pick side that is shorter (in other words, swap if 0 is longer)
            if (lengths[0] > lengths[1])
            {
                Vector3 temp1 = bendPoints[0];
                bendPoints[0] = bendPoints[1];
                bendPoints[1] = temp1;
                float temp2 = lengths[0];
                lengths[0] = lengths[1];
                lengths[1] = temp2;
                List<float> temp3 = angles[0];
                angles[0] = angles[1];
                angles[1] = temp3;
                List<MinReachRegion> temp4 = regions[0];
                regions[0] = regions[1];
                regions[1] = temp4;
                shorterSide = IntervalSide.Max;
            }
            IntervalSide longerSide = (shorterSide == IntervalSide.Min) ? IntervalSide.Max : IntervalSide.Min;

            // chain with two parts
            if (range.Size == 2)
            {
                this.angles = new float[3]
                {
                    IK2.Angle(ik.FollowChain(range, shorterSide, out _)),
                    IK2.Angle(ik.FollowChain(new IndexRangeK(range.startN, range.startN + 1, range.endN), longerSide, shorterSide)),
                    IK2.Angle(ik.FollowChain(range, longerSide, out _))
                };
                IntervalSide bendSide = (shorterSide == 0) ? IntervalSide.Max : IntervalSide.Min;
                this.regions = new MinReachRegion[2] { new MinReachRegion(bendSide, range.startN), new MinReachRegion(bendSide, range.startN + 1) };
                return;
            }

            // -------------------------
            // chain with multiple parts

            Vector3[] prevBends = new Vector3[2] { bendPoints[0], bendPoints[1] };
            Vector3[] lastBends = new Vector3[2] { bendPoints[0], bendPoints[1] };
            int[] k = new int[2] { range.startN + 1, range.startN + 1 };

            // start bending to other side starting with shorter side
            int activeIndex = 0;
            int inactiveIndex = 1;
            IntervalSide activeSide = shorterSide;
            IntervalSide inactiveSide = longerSide;
            // first part works
            while (true)
            {
                bool exit = false;
                // radius of above interval at end point that is more bent
                Interval activeInterval;
                Interval inactiveInterval;
                do
                {
                    // calculate the bend position (Forward Kinematics)
                    Vector2 reverseBend = ik.FollowChain(new IndexRangeK(range.startN, k[activeIndex], range.endN), inactiveSide, activeSide);
                    Vector2 debugCenter = ik.FollowChain(new IndexRange(0, k[activeIndex] - 1), inactiveSide, out _);
                    float debugRadius = ik.FollowChain(new IndexRange(k[activeIndex] - 1, range.endN), activeSide, out _).magnitude;
                    // update bend history
                    prevBends[activeIndex] = lastBends[activeIndex];
                    lastBends[activeIndex] = reverseBend;
                    // arc angle interval from prevBend to lastBend
                    activeInterval = IK2.AngleInterval(prevBends[activeIndex], lastBends[activeIndex]);
                    inactiveInterval = IK2.AngleInterval(prevBends[inactiveIndex], lastBends[inactiveIndex]);

                    angles[activeIndex].Insert(insertI[activeIndex], activeInterval[(int)activeSide]);
                    regions[activeIndex].Insert(insertI[activeIndex], new MinReachRegion(inactiveSide, k[activeIndex] - 1));
                    insertI[activeIndex]++;
                    k[activeIndex]++;

                    // infinite loop protection while still writing this function
                    if (k[0] > 100 || k[1] > 100)
                    {
                        exit = true;
                        break;
                    }
                    Circle otherCircle = new Circle(
                        ik.FollowChain(new IndexRange(range.startN, k[inactiveIndex] - 1), activeSide, out _),
                        ik.FollowChain(new IndexRange(k[inactiveIndex] - 1, range.endN), inactiveSide, out _).magnitude);
                    if (Vector2.Distance(otherCircle.center, reverseBend) > otherCircle.radius)
                    {
                        break;
                    }
                    if (k[activeIndex] > range.endN)
                    {
                        angles[activeIndex].Insert(insertI[activeIndex], activeInterval[(int)inactiveSide]);
                        k[activeIndex]++;
                        break;
                    }
                } while (true); // not if at end activeLength < inactiveLength, but if at any time it crossed over inactiveLength

                if (k[activeIndex] > range.endN)
                {
                    break;
                }

                activeIndex = inactiveIndex;
                inactiveIndex = (inactiveIndex + 1) % 2;
                activeSide = (activeIndex == 0) ? shorterSide : longerSide;
                inactiveSide = (activeIndex == 0) ? longerSide : shorterSide;
                // bend only the new side until there is intersection, otherwise continue with normal bends on new side
                if (Interval.AreOverlapping(activeInterval, inactiveInterval))
                {
                    break;
                }

                if (exit)
                {
                    break;
                }
            }

            if (k[activeIndex] <= range.endN)
            {
                // bend last side until there is an intersection
                int loopGuard = 100;
                do
                {
                    // continue on other side because max bend has been reached
                    if (k[activeIndex] > range.endN || k[inactiveIndex] > range.endN + 1)
                    {
                        activeIndex = inactiveIndex;
                        inactiveIndex = (inactiveIndex + 1) % 2;
                        k[activeIndex] = Mathf.Min(k[activeIndex], range.endN);
                        k[inactiveIndex] = Mathf.Min(k[inactiveIndex], range.endN + 1);
                    }
                    // centers and radiuses of last bends circles on both sides
                    Vector2 c1 = ik.FollowChain(new IndexRange(range.startN, k[activeIndex] - 2), inactiveSide, out _);
                    Vector2 c2 = ik.FollowChain(new IndexRange(range.startN, k[inactiveIndex] - 2), activeSide, out _);
                    Vector2 r1 = ik.FollowChain(new IndexRange(k[activeIndex] - 2, range.endN), activeSide, out _);
                    Vector2 r2 = ik.FollowChain(new IndexRange(k[inactiveIndex] - 2, range.endN), inactiveSide, out _);
                    
                    CircleIntersection intersection = new Circle(c1, r1.magnitude).IntersectWithCircle(new Circle(c2, r2.magnitude));

                    // there is an intersection in valid angle intervals
                    Interval arcIntersection = Interval.Overlap(new Interval(IK2.Angle(prevBends[inactiveIndex]), IK2.Angle(lastBends[inactiveIndex])),
                        new Interval(IK2.Angle(prevBends[activeIndex]), IK2.Angle(lastBends[activeIndex])));
                    if (arcIntersection.Contains(IK2.Angle(intersection.p1)))
                    {
                        angles[activeIndex].Insert(insertI[activeIndex], IK2.Angle(intersection.p1));
                        insertI[activeIndex]++;
                        break;
                    }
                    // do a normal bend
                    else
                    {
                        Vector2 reverseBend = ik.FollowChain(new IndexRangeK(range.startN, k[activeIndex], range.endN), inactiveSide, activeSide);
                        prevBends[activeIndex] = lastBends[activeIndex];
                        lastBends[activeIndex] = reverseBend;
                        Interval activeInterval = IK2.AngleInterval(prevBends[activeIndex], lastBends[activeIndex]);

                        //Debug.DrawAngle(Vector3.zero, activeInterval[inactiveIndex], 10, Color.white, 10);
                        angles[activeIndex].Insert(insertI[activeIndex], activeInterval[inactiveIndex]);
                        regions[activeIndex].Insert(insertI[activeIndex], new MinReachRegion(inactiveSide, k[activeIndex] - 1));
                        insertI[activeIndex]++;
                        k[activeIndex]++;

                        loopGuard--;
                    }
                } while (loopGuard > 0);
            }

            // save calculated values
            this.angles = new float[angles[0].Count + angles[1].Count];
            this.regions = new MinReachRegion[regions[0].Count + regions[1].Count];

            for (int i = angles[1].Count - 1; i >= 0; i--)
            {
                this.angles[angles[0].Count + (angles[1].Count - 1 - i)] = angles[1][i];
            }
            if (shorterSide == IntervalSide.Min)
            {
                for (int i = 0; i < angles[0].Count; i++)
                {
                    this.angles[i] = angles[0][i];
                }
                for (int i = angles[1].Count - 1; i >= 0; i--)
                {
                    this.angles[angles[0].Count + (angles[1].Count - 1 - i)] = angles[1][i];
                }
                for (int i = 0; i < regions[0].Count; i++)
                {
                    this.regions[i] = regions[0][i];
                }
                for (int i = regions[1].Count - 1; i >= 0; i--)
                {
                    this.regions[regions[0].Count + (regions[1].Count - 1 - i)] = regions[1][i];
                }
            }
            else
            {
                for (int i = 0; i < angles[0].Count; i++)
                {
                    this.angles[i] = angles[0][i];
                }
                for (int i = angles[1].Count - 1; i >= 0; i--)
                {
                    this.angles[angles[0].Count + (angles[1].Count - 1 - i)] = angles[1][i];
                }
                for (int i = 0; i < regions[0].Count; i++)
                {
                    this.regions[i] = regions[0][i];
                }
                for (int i = regions[1].Count - 1; i >= 0; i--)
                {
                    this.regions[regions[0].Count + (regions[1].Count - 1 - i)] = regions[1][i];
                }
            }
        }

        public float Get(IK2 ik, float angle)
        {
            for (int i = 0; i < regions.Length; i++)
            {
                if (new Interval(angles[i], angles[i + 1]).Contains(angle))
                {
                    float angleSum = 0;
                    Vector2 center = Vector2.zero;
                    for (int j = range.startN; j < regions[i].MoveChainIndex; j++)
                    {
                        angleSum += ik[j].angleLimit[(int)regions[i].StartSide];
                        center += IK2.Vector(angleSum, ik[j].length);
                    }
                    angleSum = 0;
                    Vector2 movingPart = Vector2.zero;
                    for (int j = regions[i].MoveChainIndex; j < range.endN; j++)
                    {
                        angleSum += ik[j].angleLimit[(int)regions[i].ContinueSide];
                        movingPart += IK2.Vector(angleSum, ik[j].length);
                    }

                    return new Circle(center, movingPart.magnitude).IntersectWithLine(angle).p1.magnitude; // TODO: why is everywhere p1 used???
                }
            }
            return float.NaN;
        }

        public Circle GetCircle(IK2 ik, int i)
        {
            Vector2 center = ik.FollowChain(new IndexRange(range.startN, regions[i].MoveChainIndex), regions[i].StartSide, out _);
            float radius = ik.FollowChain(new IndexRange(regions[i].MoveChainIndex, range.endN), regions[i].ContinueSide, out _).magnitude;
            return new Circle(center, radius);
        }
        public Interval GetCircleInterval(IK2 ik, int i)
        {
            Vector2 center = ik.FollowChain(new IndexRange(range.startN, regions[i].MoveChainIndex), regions[i].StartSide, out _);
            Vector2 v1 = IK2.Vector(angles[i], Get(ik, angles[i]));
            Vector2 v2 = IK2.Vector(angles[i + 1], Get(ik, angles[i + 1]));
            return new Interval(IK2.Angle(v1 - center), IK2.Angle(v2 - center));
        }

        public int AngleCount { get { return angles.Length; } }
        public Interval ValidInterval { get { return new Interval(angles[0], angles[angles.Length - 1]); } }
    }
    /// <summary>
    /// Precalculated data from which max reach at an angle is calculated
    /// </summary>
    [System.Serializable]
    public class MaxReachData
    {
        /// <summary>
        /// At index i is saved the angle interval between full bends to both sides using first i+1 segments
        /// </summary>
        private Interval[] angleSums;
        /// <summary>
        /// For which segments this data is relevant
        /// </summary>
        private IndexRange range;

        public MaxReachData(IK2 ik, IndexRange range)
        {
            this.range = new IndexRange(range.startN, range.endN);

            angleSums = new Interval[range.Size];
            for (int i = range.startN; i < range.endN; i++)
            {
                angleSums[i - range.startN] = IK2.AngleInterval(
                    ik.FollowChain(new IndexRangeK(range.startN, i + 1, range.endN), IntervalSide.Min, IntervalSide.Center),
                    ik.FollowChain(new IndexRangeK(range.startN, i + 1, range.endN), IntervalSide.Max, IntervalSide.Center));
            }
        }

        /// <summary>
        /// Get max reach of IK chain at an angle.
        /// <para>Wrapper for calling recursive version <see cref="Get(IK2, float, int, float[])"/></para>
        /// </summary>
        /// <param name="ik">Reference to IK chain</param>
        /// <param name="angle">Angle to sample max reach at</param>
        /// <returns>Max reach at angle</returns>
        public float Get(IK2 ik, float angle)
        {
            float[] lengths = new float[range.Size];
            for (int i = range.startN; i < range.endN; i++)
            {
                lengths[i - range.startN] = ik[i].length;
            }
            return Get(ik, angle, range.endN, lengths);
        }
        /// <summary>
        /// Get max reach of IK chain at an angle for segments inside <see cref="range"/> by bending chains until <see cref="k"/> (excluded) and leaving the rest extended
        /// <para>In recursive call lower k, in other words: first check longer chains and go toward shorter ones</para>
        /// </summary>
        /// <param name="ik">Reference to IK chain</param>
        /// <param name="angle">Angle to sample max reach at</param>
        /// <param name="k">Use segments inside <see cref="range"/> until <see cref="k"/> (excluded)</param>
        /// <param name="lengths">Lengths of the used segments</param>
        /// <returns></returns>
        private float Get(IK2 ik, float angle, int k, float[] lengths)
        {
            // to prevent infinite loops (just in case though it shouldn't happen)
            if (k <= range.startN) { return 0; }
            // single chain part
            else if (k - range.startN == 1) { return lengths[0]; }

            if (angleSums[k - 2 - range.startN].Contains(angle))
            {
                float[] newLengths = new float[lengths.Length - 1];
                for (int i = 0; i < newLengths.Length; i++)
                {
                    newLengths[i] = lengths[i];
                }
                newLengths[newLengths.Length - 1] += lengths[newLengths.Length];
                return Get(ik, angle, k - 1, newLengths);
            }
            else if (new Interval(angleSums[k - 1 - range.startN].min, angleSums[k - 2 - range.startN].min).Contains(angle))
            {
                Circle circle = new Circle(ik.FollowChain(new IndexRange(range.startN, k - 1), IntervalSide.Min, out _), lengths[k - 1 - range.startN]);
                return circle.IntersectWithLine(angle).p1.magnitude; // TODO: why is always p1 picked???
            }
            else if (new Interval(angleSums[k - 2 - range.startN].max, angleSums[k - 1 - range.startN].max).Contains(angle))
            {
                Circle circle = new Circle(ik.FollowChain(new IndexRange(range.startN, k - 1), IntervalSide.Max, out _), lengths[k - 1 - range.startN]);
                return circle.IntersectWithLine(angle).p1.magnitude; // TODO: why is always p1 picked???
            }
            else
            {
                return float.NaN;
            }
        }

        public Circle GetCircle(IK2 ik, int i)
        {
            int index = range.startN + Mathf.Abs(i); // from ... -2, -1, 0, 1, 2, ... to actual segment index
            Vector2 center = ik.FollowChain(new IndexRange(range.startN, index), (i < 0) ? IntervalSide.Min : IntervalSide.Max, out _);
            float radius = 0;
            for (int k = index; k < range.endN; k++)
            {
                radius += ik[k].length;
            }
            return new Circle(center, radius);
        }
        public Interval GetCircleInterval(IK2 ik, int i)
        {
            if (i == 0)
            {
                return angleSums[0];
            }
            else
            {
                IntervalSide side = (i < 0) ? IntervalSide.Min : IntervalSide.Max;
                int index = range.startN + Mathf.Abs(i); // from ... -2, -1, 0, 1, 2, ... to actual segment index
                Vector2 center = ik.FollowChain(new IndexRange(range.startN, index), side, out _);
                Vector2 v1 = ik.FollowChain(new IndexRangeK(range.startN, index, range.endN), side, IntervalSide.Center);
                Vector2 v2 = ik.FollowChain(new IndexRangeK(range.startN, index + 1, range.endN), side, IntervalSide.Center);
                v1 = v1 - center;
                v2 = v2 - center;
                return IK2.AngleInterval(v1, v2);
            }
        }

        public int AngleCount { get { return angleSums.Length; } }
        public int CircleCount { get { return angleSums.Length * 2 - 1; } }
        public Interval ValidInterval { get { return angleSums[angleSums.Length - 1]; } }
    }
    [System.Serializable]
    public struct MinReachRegion
    {
        private IntervalSide startSide;
        private int moveChainIndex;

        public MinReachRegion(IntervalSide startSide, int moveChainIndex)
        {
            this.startSide = startSide;
            this.moveChainIndex = moveChainIndex;
        }

        public override string ToString()
        {
            return $"(Start side = {StartSide}, continue side = {ContinueSide}, move chain index = {moveChainIndex})";
        }

        public IntervalSide StartSide { get { return startSide; } }
        public IntervalSide ContinueSide { get { return (startSide == IntervalSide.Min) ? IntervalSide.Max : IntervalSide.Min; } }
        public int MoveChainIndex { get { return moveChainIndex; } }
    }
}
