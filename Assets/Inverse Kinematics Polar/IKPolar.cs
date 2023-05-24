using UnityEngine;
using System.Collections.Generic;
using Unity.VisualScripting;
using static UnityEditor.Progress;
#if UNITY_EDITOR
using UnityEditor;
#endif

namespace InverseKinematicsPolarTest
{
    public enum IntervalSide : byte { Min = 0, Max = 1, Center = 2 }

    public static class IKUtility
    {
        public static Vector2 ToVector(float angle, float radius)
        {
            return new Vector2(Mathf.Cos(angle * Mathf.Deg2Rad), Mathf.Sin(angle * Mathf.Deg2Rad)) * radius;
        }
        public static float ToAngle(Vector2 point)
        {
            return Vector2.SignedAngle(Vector2.right, point);
        }
        public static Interval AngleInterval(Vector2 p1, Vector2 p2)
        {
            return new Interval(ToAngle(p1), ToAngle(p2));
        }
        public static Vector2 FollowChain(IKSegment[] segments, IndexRangeK fromToSplit, IntervalSide limitIndex1, IntervalSide limitIndex2)
        {
            Vector2 p1 = FollowChain(segments, new IndexRange(fromToSplit.startN, fromToSplit.k), limitIndex1, out float endAngle1);
            Vector2 p2 = FollowChain(segments, new IndexRange(fromToSplit.k, fromToSplit.endN), limitIndex2);
            return p1 + (Vector2)(Quaternion.Euler(0, 0, endAngle1) * p2);
        }
        public static Vector2 FollowChain(IKSegment[] segments, IndexRange fromTo, IntervalSide limitIndex)
        {
            return FollowChain(segments, fromTo, limitIndex, out float endAngle);
        }
        public static Vector2 FollowChain(IKSegment[] segments, IndexRange fromTo, IntervalSide limitIndex, out float endAngle)
        {
            endAngle = 0;
            Vector2 p = Vector2.zero;
            for (int i = fromTo.startN; i < fromTo.endN; i++)
            {
                if (limitIndex == IntervalSide.Min || limitIndex == IntervalSide.Max)
                {
                    endAngle += segments[i].angleLimit[(int)limitIndex];
                }
                p += ToVector(endAngle, segments[i].length);
            }
            return p;
        }

        public static void DebugFollowChain(IKSegment[] segments, IndexRangeK fromToSplit, IntervalSide limitIndex1, IntervalSide limitIndex2, float startAngle, Color color)
        {
            Vector2 p1 = FollowChain(segments, new IndexRange(fromToSplit.startN, fromToSplit.k), limitIndex1, out float endAngle);
            DebugFollowChain(segments, new IndexRange(fromToSplit.startN, fromToSplit.k), limitIndex1, startAngle, Vector2.zero, color, out float endAngle1);
            DebugFollowChain(segments, new IndexRange(fromToSplit.k, fromToSplit.endN), limitIndex2, endAngle1, p1, color);
        }
        public static void DebugFollowChain(IKSegment[] segments, IndexRange fromTo, IntervalSide limitIndex, float startAngle, Vector2 origin, Color color)
        {
            DebugFollowChain(segments, fromTo, limitIndex, startAngle, origin, color, out float endAngle);
        }
        public static void DebugFollowChain(IKSegment[] segments, IndexRange fromTo, IntervalSide limitIndex, float startAngle, Vector2 origin, Color color, out float endAngle)
        {
            endAngle = startAngle;
            Vector2 p = origin;
            for (int i = fromTo.startN; i < fromTo.endN; i++)
            {
                if (limitIndex == IntervalSide.Min || limitIndex == IntervalSide.Max)
                {
                    endAngle += segments[i].angleLimit[(int)limitIndex];
                }
                Vector2 nextChain = ToVector(endAngle, segments[i].length);
                Debug.DrawLine(p, p + nextChain, color, 10);
                p += nextChain;
            }
        }

        public static CircleIntersection TwoCircleIntersection(Circle c1, Circle c2)
        {
            return TwoCircleIntersection(c1.center, c2.center, c1.radius, c2.radius);
        }
        public static CircleIntersection TwoCircleIntersection(Circle c1, Vector2 c2, float r2)
        {
            return TwoCircleIntersection(c1.center, c2, c1.radius, r2);
        }
        public static CircleIntersection TwoCircleIntersection(Vector2 c1, Vector2 c2, float r1, float r2)
        {
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
                CircleIntersection intersection = IntersectionLineCircle(90, new Vector2(x - c1.x, 0), r1);
                intersection.p1 = new Vector2(x, intersection.p1.y + c1.y);
                intersection.p2 = new Vector2(x, intersection.p2.y + c1.y);
                return intersection;
                //return new CircleIntersection(new Vector2(x, c1.y), new Vector2(float.NaN, float.NaN));
            }
            else
            {
                return (r1 == r2) ? new CircleIntersection(c1.x, c1.y, r1, float.PositiveInfinity) : CircleIntersection.Miss;
            }
        }
        public static CircleIntersection IntersectionLineCircle(float angle, Vector2 center, float radius)
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
    }

    public class IKPolar : MonoBehaviour
    {
        public bool debugChainBend = true;
        public bool debugTargetLine = true;
        public int debugN = 4;
        [Range(0, 5)] public int debugLevel;
        public Transform target;

        public IKSegment[] segments;
        private Reach[] dataFront;
        private Reach[] dataBack;

        private Vector2 mousePos = Vector2.zero;
        private Vector2 reach = Vector2.zero;

        private void Awake()
        {
            DebugChainBend(Vector2.zero, 0, debugN, Color.green, 10);
            Precalculate();
        }

        private void Update()
        {
            for (int i = 0; i < segments.Length; i++)
            {
                if (segments[i].animateTargetAngle)
                {
                    segments[i].targetAngle = (Mathf.Sin(Time.time) + 1) / 2;
                }
            }

            if (debugChainBend)
            {
                DebugChainBend(Vector2.zero, 0, 0, Color.green, 0.01f);
            }

            if (target == null)
            {
                mousePos = Camera.main.ScreenToWorldPoint(new Vector3(Input.mousePosition.x, Input.mousePosition.y, -Camera.main.transform.position.z));
            }
            else
            {
                mousePos = target.position;
            }
            mousePos = ClampTarget(mousePos);
            float angle = Mathf.Asin(mousePos.y / mousePos.magnitude) * Mathf.Rad2Deg;


            Interval reachInterval = dataFront[0].Get(segments, angle);
            reach = new Vector2(reachInterval.min, reachInterval.max);
            Debug.DrawRay(mousePos.normalized * reach.x, mousePos.normalized * Mathf.Max(reach.y - reach.x, Mathf.Epsilon), Color.red, 0.01f);
            if (debugTargetLine)
            {
                Debug.DrawRay(Vector3.zero, mousePos, Color.yellow, 0.01f);
            }

            Debug.DrawReach(dataFront[0], segments, Color.black);
            Align2(Vector3.zero, mousePos, -45, segments.Length);
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

        private void DebugChainBend(Vector2 pos, float angle, int startN, Color color, float duration)
        {
            DebugChainBend(pos, Quaternion.Euler(0, 0, angle), startN, color, duration);
        }
        private void DebugChainBend(Vector2 pos, Quaternion rot, int startN, Color color, float duration)
        {
            // draw min/max bend
            Vector2 pMin = pos;
            Vector2 pMax = pos;

            Interval angleSum = new Interval(0, 0);
            for (int i = startN; i < segments.Length; i++)
            {
                angleSum.min += segments[i].angleLimit[0];
                angleSum.max += segments[i].angleLimit[1];
                Vector2 dirMin = rot * new Vector2(Mathf.Cos(angleSum.min * Mathf.Deg2Rad), Mathf.Sin(angleSum.min * Mathf.Deg2Rad)) * segments[i].length;
                Vector2 dirMax = rot * new Vector2(Mathf.Cos(angleSum.max * Mathf.Deg2Rad), Mathf.Sin(angleSum.max * Mathf.Deg2Rad)) * segments[i].length;
                Debug.DrawRay(pMin, dirMin, color, duration);
                Debug.DrawRay(pMax, dirMax, color, duration);
                pMin += dirMin;
                pMax += dirMax;
            }
        }

        private void Precalculate()
        {
            dataFront = new Reach[segments.Length];
            int i = 0;
            for (int k = segments.Length; k > 0; k--)
            {
                dataFront[i] = new Reach(segments, 0, k);
                i++;
            }
            dataBack = new Reach[segments.Length];
            for (int j = 0; j < segments.Length; j++)
            {
                dataBack[j] = new Reach(segments, j, segments.Length);
            }
        }

        private Vector2 ClampTarget(Vector2 target)
        {
            float angle = IKUtility.ToAngle(target);
            float radius = target.magnitude;

            angle = Interval.ClampTo(angle, dataFront[0].ValidInterval);
            radius = Interval.ClampTo(radius, dataFront[0].Get(segments, angle));
            return radius * new Vector2(Mathf.Cos(angle * Mathf.Deg2Rad), Mathf.Sin(angle * Mathf.Deg2Rad));
        }
        public void Align(Vector2 target, float angle, int n)
        {
            if (n == 1)
            {
                Debug.DrawLine(transform.position, target, Color.cyan);
                return;
            }
            else if (n == 2)
            {
                float a = segments[0].length;
                float b = target.magnitude;
                float c = segments[1].length;
                float addAngle = Mathf.Acos(Mathf.Clamp((a * a + b * b - c * c) / (2 * a * b), -1, 1)) * Mathf.Rad2Deg;
                Vector2 p;
                Vector2 dir;
                p = IKUtility.ToVector(IKUtility.ToAngle(target) + addAngle, segments[0].length);
                if (segments[0].angleLimit.Contains(IKUtility.ToAngle(p)))
                {
                    dir = (target - p).normalized * segments[1].length;
                    float angleDif = angle - IKUtility.ToAngle(dir);
                    if (segments[2].angleLimit.Contains(angleDif))
                    {
                        Debug.DrawLine(Vector2.zero, p, Color.yellow);
                        Debug.DrawLine(p, p + dir, Color.yellow);
                    }
                    else
                    {
                        Debug.DrawLine(Vector2.zero, p, Color.red);
                        Debug.DrawLine(p, p + dir, Color.red);
                    }
                }

                p = IKUtility.ToVector(IKUtility.ToAngle(target) - addAngle, segments[0].length);
                if (segments[0].angleLimit.Contains(IKUtility.ToAngle(p)))
                {
                    dir = (target - p).normalized * segments[1].length;
                    float angleDif = angle - IKUtility.ToAngle(dir);
                    if (segments[2].angleLimit.Contains(angleDif))
                    {
                        Debug.DrawLine(Vector2.zero, p, Color.blue);
                        Debug.DrawLine(p, p + dir, Color.blue);
                    }
                    else
                    {
                        Debug.DrawLine(Vector2.zero, p, Color.red);
                        Debug.DrawLine(p, p + dir, Color.red);
                    }
                }

                return;
            }

            Circle targetCircle = new Circle(target, segments[n - 1].length);

            Debug.DrawCircle(targetCircle, 64, Color.magenta);
            Debug.DrawLine(target, target + IKUtility.ToVector(angle, segments[n - 1].length) / 2, Color.magenta);

            Reach reach = GetData(n - 1);
            Debug.DrawReach(reach, segments, Color.green);

            List<Vector2> minPoints = new List<Vector2>();
            List<Vector2> maxPoints = new List<Vector2>();
            foreach (var circleInterval in reach.GetMinCircles(segments))
            {
                Circle circle = circleInterval.Key;
                CircleIntersection Inter = IKUtility.TwoCircleIntersection(circle, targetCircle);

                switch (Inter.Type)
                {
                    case CircleIntersection.Variant.Touching:
                        minPoints.Add(Inter.p1);
                        break;
                    case CircleIntersection.Variant.Intersect:
                        float a1 = IKUtility.ToAngle(Inter.p1 - circleInterval.Key.center);
                        float a2 = IKUtility.ToAngle(Inter.p2 - circleInterval.Key.center);
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
            foreach (var circleInterval in reach.GetMaxCircles(segments))
            {
                Circle circle = circleInterval.Key;
                CircleIntersection Inter = IKUtility.TwoCircleIntersection(circle, targetCircle);

                switch (Inter.Type)
                {
                    case CircleIntersection.Variant.Touching:
                        maxPoints.Add(Inter.p1);
                        break;
                    case CircleIntersection.Variant.Intersect:
                        float a1 = IKUtility.ToAngle(Inter.p1 - circleInterval.Key.center);
                        float a2 = IKUtility.ToAngle(Inter.p2 - circleInterval.Key.center);
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

            if (maxPoints.Count == 0)
            {
                float a = IKUtility.ToAngle(target);
                a = Interval.ClampTo(a, reach.ValidInterval);
                maxPoints.Add(IKUtility.ToVector(a, reach.Get(segments, a).max));
            }
            else if (maxPoints.Count == 1)
            {
                Vector2 v = IKUtility.FollowChain(segments, new IndexRange(0, n - 1), IntervalSide.Min);
                if (Vector2.Distance(targetCircle.center, v) <= targetCircle.radius)
                {
                    maxPoints.Add(v);
                }
                v = IKUtility.FollowChain(segments, new IndexRange(0, n - 1), IntervalSide.Max);
                if (Vector2.Distance(targetCircle.center, v) <= targetCircle.radius)
                {
                    maxPoints.Add(v);
                }
            }
            else
            {
                if (IKUtility.ToAngle(maxPoints[0]) > IKUtility.ToAngle(maxPoints[1]))
                {
                    Vector2 temp = maxPoints[0];
                    maxPoints[0] = maxPoints[1];
                    maxPoints[1] = temp;
                }
            }
            minPoints.Sort((a, b) => IKUtility.ToAngle(a).CompareTo(IKUtility.ToAngle(b)));

            Interval newTargetInterval = new Interval(IKUtility.ToAngle(maxPoints[0]), IKUtility.ToAngle(maxPoints[1 % maxPoints.Count]));

            List<Interval> intervals = new List<Interval>();
            if (maxPoints.Count == 1)
            {
                intervals.Add(new Interval(IKUtility.ToAngle(maxPoints[0]), IKUtility.ToAngle(maxPoints[0])));
            }
            else
            {
                if (minPoints.Count == 0)
                {
                    intervals.Add(new Interval(IKUtility.ToAngle(maxPoints[0]), IKUtility.ToAngle(maxPoints[1])));
                }
                else if (minPoints.Count % 2 == 0)
                {
                    intervals.Add(new Interval(IKUtility.ToAngle(maxPoints[0]), IKUtility.ToAngle(minPoints[0])));
                    for (int i = 1; i < minPoints.Count - 1; i += 2)
                    {
                        intervals.Add(new Interval(IKUtility.ToAngle(minPoints[i]), IKUtility.ToAngle(minPoints[i + 1])));
                    }
                    intervals.Add(new Interval(IKUtility.ToAngle(minPoints[minPoints.Count - 1]), IKUtility.ToAngle(maxPoints[1])));
                }
                else
                {
                    if (maxPoints[0].magnitude < maxPoints[1].magnitude)
                    {
                        for (int i = 0; i < minPoints.Count - 1; i += 2)
                        {
                            intervals.Add(new Interval(IKUtility.ToAngle(minPoints[i]), IKUtility.ToAngle(minPoints[i + 1])));
                        }
                        intervals.Add(new Interval(IKUtility.ToAngle(minPoints[minPoints.Count - 1]), IKUtility.ToAngle(maxPoints[1])));
                    }
                    else
                    {
                        intervals.Add(new Interval(IKUtility.ToAngle(maxPoints[0]), IKUtility.ToAngle(minPoints[0])));
                        for (int i = 0; i < minPoints.Count - 1; i += 2)
                        {
                            intervals.Add(new Interval(IKUtility.ToAngle(minPoints[i]), IKUtility.ToAngle(minPoints[i + 1])));
                        }
                    }
                }
            }
            MultiInterval mi = new MultiInterval(intervals);
            float newAngle = mi.GetConnected(segments[n - 2].targetAngle);
            CircleIntersection inter = IKUtility.IntersectionLineCircle(newAngle, targetCircle.center, targetCircle.radius);
            if (inter.Type == CircleIntersection.Variant.Miss)
            {
                Debug.DrawCircle(targetCircle, 64, Color.cyan);
                Debug.DrawLine(Vector3.zero, IKUtility.ToVector(newAngle, 10), Color.cyan);
            }
            Vector2 newTarget = IKUtility.ToVector(newAngle, Mathf.Min(inter.p1.magnitude, inter.p2.magnitude));

            Debug.DrawLine(targetCircle.center, newTarget, Color.yellow);

            Align(newTarget, IKUtility.ToAngle(targetCircle.center - newTarget), n - 1);
        }


        public float debug = 0;


        public void Align2(Vector3 origin, Vector2 target, float angle, int n, float prevAngle = 0)
        {
            if (n == 1)
            {
                Debug.DrawRay(origin, (target - (Vector2)origin).normalized * segments[segments.Length - 1].length, Color.cyan);
                return;
            }
            /*else if (n == 2)
            {
                return;
                float a = segments[0].length;
                float b = target.magnitude;
                float c = segments[1].length;
                float addAngle = Mathf.Acos(Mathf.Clamp((a * a + b * b - c * c) / (2 * a * b), -1, 1)) * Mathf.Rad2Deg;
                Vector2 p;
                Vector2 dir;
                p = IKUtility.ToVector(IKUtility.ToAngle(target) + addAngle, segments[0].length);
                if (segments[0].angleLimit.Contains(IKUtility.ToAngle(p)))
                {
                    dir = (target - p).normalized * segments[1].length;
                    float angleDif = angle - IKUtility.ToAngle(dir);
                    if (segments[2].angleLimit.Contains(angleDif))
                    {
                        Debug.DrawLine(Vector2.zero, p, Color.cyan);
                        Debug.DrawLine(p, p + dir, Color.cyan);
                    }
                    else
                    {
                        Debug.DrawLine(Vector2.zero, p, Color.green);
                        Debug.DrawLine(p, p + dir, Color.green);
                    }
                }

                p = IKUtility.ToVector(IKUtility.ToAngle(target) - addAngle, segments[0].length);
                if (segments[0].angleLimit.Contains(IKUtility.ToAngle(p)))
                {
                    dir = (target - p).normalized * segments[1].length;
                    float angleDif = angle - IKUtility.ToAngle(dir);
                    if (segments[2].angleLimit.Contains(angleDif))
                    {
                        Debug.DrawLine(Vector2.zero, p, Color.cyan);
                        Debug.DrawLine(p, p + dir, Color.cyan);
                    }
                    else
                    {
                        Debug.DrawLine(Vector2.zero, p, Color.blue);
                        Debug.DrawLine(p, p + dir, Color.blue);
                    }
                }

                return;
            }*/
            Debug.DrawPoint(target, Color.white, 3);
            bool shouldDebug = debugLevel == n;

            Debug.Log(n);
            int segmentI = segments.Length - n;
            Circle targetCircle = new Circle(new Vector2(-segments[segmentI].length, 0), Vector3.Distance(origin, target));

            Reach smallerReach = dataBack[segmentI + 1];
            float amin = segments[segmentI].angleLimit.min;
            float amax = segments[segmentI].angleLimit.max;
            Vector3 c1 = IKUtility.ToVector(amin, segments[segmentI].length);
            Vector3 c2 = IKUtility.ToVector(amax, segments[segmentI].length);
            //Debug.DrawReach(c1, smallerReach, segments, Color.cyan, amin);
            //Debug.DrawReach(c2, smallerReach, segments, Color.yellow, amax);
            float rotDebugAngle = prevAngle + ((segments[segmentI].targetAngle * 2) - 1) * 180;
            Vector2 localTarget = new Vector2(float.NaN, float.NaN);
            Vector2 toTarget = target - (Vector2)origin;
            float ta = Mathf.DeltaAngle(prevAngle, IKUtility.ToAngle(toTarget));
            localTarget = IKUtility.ToVector(ta, toTarget.magnitude);
            if (shouldDebug)
            {
                Debug.DrawCircle(new Circle(targetCircle.center + new Vector2(segments[segmentI].length, 0), targetCircle.radius), 64, Color.black);
                /*Debug.DrawReach(origin + Quaternion.Euler(0, 0, rotDebugAngle) * new Vector3(segments[segmentI].length, 0), 
                    smallerReach, segments, Color.red, rotDebugAngle);*/

                //Debug.DrawLine(origin, target, Color.white);
                Debug.DrawPoint(targetCircle.center + localTarget, new Color(1.0f, 0.75f, 0.0f));

                Debug.DrawCircle(targetCircle, 64, Color.red);
                Debug.DrawReach(targetCircle.center + IKUtility.ToVector(prevAngle, segments[segmentI].length), smallerReach, segments, Color.yellow, prevAngle);
                Debug.DrawReach(smallerReach, segments, Color.magenta);
            }

            // get intersection intervals between target circle and reach (facing directly right)
            List<Vector2> minPoints = new List<Vector2>();
            List<Vector2> maxPoints = new List<Vector2>();
            foreach (var circleInterval in smallerReach.GetMinCircles(segments))
            {
                Circle circle = circleInterval.Key;
                CircleIntersection Inter = IKUtility.TwoCircleIntersection(circle, targetCircle);

                switch (Inter.Type)
                {
                    case CircleIntersection.Variant.Touching:
                        minPoints.Add(Inter.p1 + new Vector2(segments[segmentI].length, 0));
                        break;
                    case CircleIntersection.Variant.Intersect:
                        float a1 = IKUtility.ToAngle(Inter.p1 - circle.center);
                        float a2 = IKUtility.ToAngle(Inter.p2 - circle.center);
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
            foreach (var circleInterval in smallerReach.GetMaxCircles(segments))
            {
                Circle circle = circleInterval.Key;
                CircleIntersection Inter = IKUtility.TwoCircleIntersection(circle, targetCircle);

                switch (Inter.Type)
                {
                    case CircleIntersection.Variant.Touching:
                        maxPoints.Add(Inter.p1 + new Vector2(segments[segmentI].length, 0));
                        break;
                    case CircleIntersection.Variant.Intersect:
                        float a1 = IKUtility.ToAngle(Inter.p1 - circle.center);
                        float a2 = IKUtility.ToAngle(Inter.p2 - circle.center);
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

            if (shouldDebug)
            {
                Debug.Log(minPoints.Count + " .. " + maxPoints.Count);
                foreach (var item in minPoints)
                {
                    Debug.DrawPoint(item, Color.red);
                }
                foreach (var item in maxPoints)
                {
                    Debug.DrawPoint(item, Color.green);
                }
            }

            #region multiinterval magic

            if (maxPoints.Count == 0)
            {
                float a = IKUtility.ToAngle(target);
                a = Interval.ClampTo(a, smallerReach.ValidInterval);
                maxPoints.Add(IKUtility.ToVector(a, smallerReach.Get(segments, a).max));
            }
            else if (maxPoints.Count == 1)
            {
                Vector2 v = IKUtility.FollowChain(segments, new IndexRange(0, n - 1), IntervalSide.Min);
                if (Vector2.Distance(targetCircle.center, v) <= targetCircle.radius)
                {
                    maxPoints.Add(v);
                }
                v = IKUtility.FollowChain(segments, new IndexRange(0, n - 1), IntervalSide.Max);
                if (Vector2.Distance(targetCircle.center, v) <= targetCircle.radius)
                {
                    maxPoints.Add(v);
                }
            }
            else
            {
                if (IKUtility.ToAngle(maxPoints[0]) > IKUtility.ToAngle(maxPoints[1]))
                {
                    Vector2 temp = maxPoints[0];
                    maxPoints[0] = maxPoints[1];
                    maxPoints[1] = temp;
                }
            }
            minPoints.Sort((a, b) => IKUtility.ToAngle(a).CompareTo(IKUtility.ToAngle(b)));

            Interval newTargetInterval = new Interval(IKUtility.ToAngle(maxPoints[0]), IKUtility.ToAngle(maxPoints[1 % maxPoints.Count]));

            List<Interval> intervals = new List<Interval>();
            if (maxPoints.Count == 1)
            {
                if (minPoints.Count == 0)
                {
                    intervals.Add(new Interval(IKUtility.ToAngle(maxPoints[0] - targetCircle.center), IKUtility.ToAngle(maxPoints[0] - targetCircle.center)));
                }
                else
                {
                    intervals.Add(new Interval(IKUtility.ToAngle(minPoints[0] - targetCircle.center), IKUtility.ToAngle(maxPoints[0] - targetCircle.center)));
                }
            }
            else
            {
                if (minPoints.Count == 0)
                {
                    intervals.Add(new Interval(IKUtility.ToAngle(maxPoints[0] - targetCircle.center), IKUtility.ToAngle(maxPoints[1] - targetCircle.center)));
                }
                else if (minPoints.Count % 2 == 0)
                {
                    intervals.Add(new Interval(IKUtility.ToAngle(maxPoints[0] - targetCircle.center), IKUtility.ToAngle(minPoints[0] - targetCircle.center)));
                    for (int i = 1; i < minPoints.Count - 1; i += 2)
                    {
                        intervals.Add(new Interval(IKUtility.ToAngle(minPoints[i] - targetCircle.center), IKUtility.ToAngle(minPoints[i + 1] - targetCircle.center)));
                    }
                    intervals.Add(new Interval(IKUtility.ToAngle(minPoints[minPoints.Count - 1] - targetCircle.center), IKUtility.ToAngle(maxPoints[1] - targetCircle.center)));
                }
                else
                {
                    if (maxPoints[0].magnitude < maxPoints[1].magnitude)
                    {
                        for (int i = 0; i < minPoints.Count - 1; i += 2)
                        {
                            intervals.Add(new Interval(IKUtility.ToAngle(minPoints[i]), IKUtility.ToAngle(minPoints[i + 1])));
                        }
                        intervals.Add(new Interval(IKUtility.ToAngle(minPoints[minPoints.Count - 1]), IKUtility.ToAngle(maxPoints[1])));
                    }
                    else
                    {
                        intervals.Add(new Interval(IKUtility.ToAngle(maxPoints[0]), IKUtility.ToAngle(minPoints[0])));
                        for (int i = 0; i < minPoints.Count - 1; i += 2)
                        {
                            intervals.Add(new Interval(IKUtility.ToAngle(minPoints[i]), IKUtility.ToAngle(minPoints[i + 1])));
                        }
                    }
                }
            }
            #endregion
            // rotate all intervals by previous part angle
            for (int i = 0; i < intervals.Count; i++)
            {
                //intervals[i] += prevAngle;
            }
            //Debug.Log(intervals[0].IsPoint());

            if (shouldDebug)
            {
                foreach (var item in intervals)
                {
                    Debug.DrawAngleInterval(targetCircle.center, item, targetCircle.radius, Color.yellow);
                }
            }

            float targetAngle;
            if (segmentI > 0)
            {
                targetAngle = IKUtility.ToAngle(target - (Vector2)(origin - (Vector3)IKUtility.ToVector(prevAngle, segments[segmentI - 1].length)));
            }
            else
            {
                targetAngle = IKUtility.ToAngle(target);
            }
            targetAngle = IKUtility.ToAngle(localTarget);
            if (shouldDebug)
            {
                //Debug.DrawAngle(targetCircle.center, targetAngle, 10, Color.green);
                Debug.DrawAngleInterval(targetCircle.center, segments[segmentI].angleLimit, 10, Color.black);
            }
            float partAngle;
            if (intervals.Count == 1 && intervals[0].IsPoint())
            {
                partAngle = prevAngle + targetAngle - intervals[0].min;
                //Debug.DrawAngle(targetCircle.center, intervals[0].min, targetCircle.radius, Color.green);
                Debug.Log(shouldDebug + "  " + n);
            }
            else
            {
                List<Interval> angleOffsets = new List<Interval>();
                for (int i = 0; i < intervals.Count; i++)
                {
                    //if (intervals[i].Contains(targetAngle)) { continue; }
                    Interval valid = Interval.Overlap(new Interval(targetAngle - intervals[i].min, targetAngle - intervals[i].max), segments[segmentI].angleLimit);
                    if (shouldDebug)
                    {
                        Debug.Log(intervals[i] + "  " + targetAngle);
                        Debug.Log(new Interval(targetAngle - intervals[i].min, targetAngle - intervals[i].max));
                        Debug.DrawAngleInterval(targetCircle.center, new Interval(targetAngle - intervals[i].min, targetAngle - intervals[i].max), 30, Color.blue);
                    }

                    if (!valid.IsEmpty())
                    {
                        angleOffsets.Add(valid);
                    }
                }
                angleOffsets.Sort((a, b) => a.min.CompareTo(b.min));

                if (shouldDebug)
                {
                    foreach (var item in angleOffsets)
                    {
                        //Debug.DrawAngleInterval(targetCircle.center, item, targetCircle.radius, Color.blue);
                    }
                }
                MultiInterval mi = new MultiInterval(angleOffsets);
                partAngle = prevAngle + mi.GetConnected(segments[segmentI].targetAngle);
            }
            if (shouldDebug)
            {
                Debug.DrawAngle(targetCircle.center, partAngle, 10, Color.green);
            }
            Vector3 partVec = IKUtility.ToVector(partAngle, segments[segmentI].length);
            Debug.DrawRay(origin, partVec, Color.cyan);
            //Debug.Log(origin + "  " + (origin + partVec) + " " + partAngle);
            /*if (n == 4)
            {
                Debug.DrawReach(origin + new Vector3(segments[segmentI].length, 0), smallerReach, segments, Color.white, partAngle);
                Debug.DrawReach(origin + partVec + new Vector3(segments[segmentI + 1].length, 0), dataBack[segmentI + 2], segments, Color.yellow, partAngle);
            }*/

            Align2(origin + partVec, target, angle, n - 1, partAngle);
        }

        private Reach GetData(int n)
        {
            return dataFront[segments.Length - n];
        }
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
    [System.Serializable]
    public class Reach
    {
        private MinReachData minData;
        private MaxReachData maxData;

        public Reach(IKSegment[] segments, int startN, int endN)
        {
            minData = new MinReachData(segments, startN, endN);
            maxData = new MaxReachData(segments, startN, endN);
        }

        public IEnumerable<KeyValuePair<Circle, Interval>> GetMinCircles(IKSegment[] segments)
        {
            for (int i = 0; i < minData.AngleCount - 1; i++)
            {
                yield return new KeyValuePair<Circle, Interval>(minData.GetCircle(segments, i), minData.GetCircleInterval(segments, i));
            }
        }
        public IEnumerable<KeyValuePair<Circle, Interval>> GetMaxCircles(IKSegment[] segments)
        {
            for (int i = -(maxData.AngleCount - 1); i < maxData.AngleCount; i++)
            {
                yield return new KeyValuePair<Circle, Interval>(maxData.GetCircle(segments, i), maxData.GetCircleInterval(segments, i));
            }
        }

        public Interval Get(IKSegment[] segments, float angle)
        {
            return new Interval(minData.Get(segments, angle), maxData.Get(segments, angle));
        }
        public Interval ValidInterval { get { return maxData.ValidInterval; } }
    }
    [System.Serializable]
    public class MinReachData
    {
        private float[] angles;
        private MinReachRegion[] regions;
        private int startN, partCount;

        public MinReachData(IKSegment[] segments, int startN, int endN)
        {
            this.startN = startN;
            this.partCount = endN;

            //Debug.Log("_ " + startN + "  " + endN);

            if (endN - startN <= 0)
            {
                return;
            }
            // only one chain part
            if (endN - startN == 1)
            {
                this.angles = new float[2] { segments[startN].angleLimit.min, segments[startN].angleLimit.max };
                this.regions = new MinReachRegion[1] { new MinReachRegion(IntervalSide.Min, startN) };
                return;
            }
            // [0]: shorter bend side
            // [1]: longer bend side
            IntervalSide shorterSide = IntervalSide.Min;
            Vector3[] bendPoints = new Vector3[2]
            {
            IKUtility.FollowChain(segments, new IndexRange(startN, endN), IntervalSide.Min),
            IKUtility.FollowChain(segments, new IndexRange(startN, endN), IntervalSide.Max)
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
            //Debug.Log("_ shorter side " + shorterSide);

            // chain with two parts
            if (endN - startN == 2)
            {
                this.angles = new float[3]
                {
                IKUtility.ToAngle(IKUtility.FollowChain(segments, new IndexRange(startN, endN), shorterSide)),
                IKUtility.ToAngle(IKUtility.FollowChain(segments, new IndexRangeK(startN, startN + 1, endN), longerSide, shorterSide)),
                IKUtility.ToAngle(IKUtility.FollowChain(segments, new IndexRange(startN, endN), longerSide))
                };
                IntervalSide bendSide = (shorterSide == 0) ? IntervalSide.Max : IntervalSide.Min;
                this.regions = new MinReachRegion[2] { new MinReachRegion(bendSide, startN), new MinReachRegion(bendSide, startN + 1) };
                return;
            }

            // -------------------------
            // chain with multiple parts

            Vector3[] prevBends = new Vector3[2] { bendPoints[0], bendPoints[1] };
            Vector3[] lastBends = new Vector3[2] { bendPoints[0], bendPoints[1] };
            int[] k = new int[2] { startN + 1, startN + 1 };

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
                    //Debug.Log(endN + "  " + k[activeIndex]);
                    //IKUtility.DebugFollowChain(segments, new IndexRangeK(startN, k[activeIndex], endN), inactiveSide, activeSide, 0, Color.yellow);
                    Vector2 reverseBend = IKUtility.FollowChain(segments, new IndexRangeK(startN, k[activeIndex], endN), inactiveSide, activeSide);
                    Vector2 debugCenter = IKUtility.FollowChain(segments, new IndexRange(0, k[activeIndex] - 1), inactiveSide);
                    float debugRadius = IKUtility.FollowChain(segments, new IndexRange(k[activeIndex] - 1, endN), activeSide).magnitude;
                    //Debug.DrawCircle(debugCenter, debugRadius, 64, Color.cyan, 10);
                    //Debug.Log(activeSide + "  " + k[activeIndex]);
                    // update bend history
                    prevBends[activeIndex] = lastBends[activeIndex];
                    lastBends[activeIndex] = reverseBend;
                    // arc angle interval from prevBend to lastBend
                    activeInterval = IKUtility.AngleInterval(prevBends[activeIndex], lastBends[activeIndex]);
                    inactiveInterval = IKUtility.AngleInterval(prevBends[inactiveIndex], lastBends[inactiveIndex]);

                    //Debug.Log("~> " + activeIndex + "  " + k[activeIndex] + " " + activeSide + "  " + activeInterval);

                    //Debug.DrawAngle(Vector3.zero, activeInterval[(int)activeSide], 10, Color.black, 10);
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
                        IKUtility.FollowChain(segments, new IndexRange(startN, k[inactiveIndex] - 1), activeSide),
                        IKUtility.FollowChain(segments, new IndexRange(k[inactiveIndex] - 1, endN), inactiveSide).magnitude);
                    //Debug.DrawCircle(otherCircle, 64, Color.cyan, 10);
                    if (Vector2.Distance(otherCircle.center, reverseBend) > otherCircle.radius)
                    {
                        break;
                    }
                    if (k[activeIndex] > endN)
                    {
                        //Debug.DrawAngle(Vector3.zero, activeInterval[(int)inactiveSide], 10, Color.black, 10);
                        angles[activeIndex].Insert(insertI[activeIndex], activeInterval[(int)inactiveSide]);
                        k[activeIndex]++;
                        break;
                    }
                } while (true); // not if at end activeLength < inactiveLength, but if at any time it crossed over inactiveLength

                if (k[activeIndex] > endN)
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

            if (k[activeIndex] <= endN)
            {
                Debug.Log("Do Part 2");
                // bend last side until there is an intersection
                int loopGuard = 100;
                do
                {
                    // continue on other side because max bend has been reached
                    if (k[activeIndex] > endN || k[inactiveIndex] > endN + 1)
                    {
                        activeIndex = inactiveIndex;
                        inactiveIndex = (inactiveIndex + 1) % 2;
                        k[activeIndex] = Mathf.Min(k[activeIndex], endN);
                        k[inactiveIndex] = Mathf.Min(k[inactiveIndex], endN + 1);
                    }
                    // centers and radiuses of last bends circles on both sides
                    Vector2 c1 = IKUtility.FollowChain(segments, new IndexRange(startN, k[activeIndex] - 2), inactiveSide);
                    Vector2 c2 = IKUtility.FollowChain(segments, new IndexRange(startN, k[inactiveIndex] - 2), activeSide);
                    Vector2 r1 = IKUtility.FollowChain(segments, new IndexRange(k[activeIndex] - 2, endN), activeSide);
                    Vector2 r2 = IKUtility.FollowChain(segments, new IndexRange(k[inactiveIndex] - 2, endN), inactiveSide);
                    //Debug.DrawCircle(c1, r1.magnitude, 128, Color.white, 10);
                    //Debug.DrawCircle(c2, r2.magnitude, 128, Color.black, 10);
                    CircleIntersection intersection = IKUtility.TwoCircleIntersection(c1, c2, r1.magnitude, r2.magnitude);

                    // there is an intersection in valid angle intervals
                    Interval arcIntersection = Interval.Overlap(new Interval(IKUtility.ToAngle(prevBends[inactiveIndex]), IKUtility.ToAngle(lastBends[inactiveIndex])),
                        new Interval(IKUtility.ToAngle(prevBends[activeIndex]), IKUtility.ToAngle(lastBends[activeIndex])));
                    //Debug.DrawAngleInterval(Vector3.zero, arcIntersection, 20, Color.magenta, 10);
                    if (arcIntersection.Contains(IKUtility.ToAngle(intersection.p1)))
                    {
                        //Debug.DrawAngle(Vector3.zero, IKUtility.ToAngle(intersection.p1), 10, Color.white, 10);
                        angles[activeIndex].Insert(insertI[activeIndex], IKUtility.ToAngle(intersection.p1));
                        insertI[activeIndex]++;
                        break;
                    }
                    // do a normal bend
                    else
                    {
                        //IKUtility.DebugFollowChain(segments, new IndexRangeK(startN, k[activeIndex], partCount), inactiveSide, activeSide, 0, new Color(1.0f, 0.5f, 0.0f));
                        Vector2 reverseBend = IKUtility.FollowChain(segments, new IndexRangeK(startN, k[activeIndex], endN), inactiveSide, activeSide);
                        prevBends[activeIndex] = lastBends[activeIndex];
                        lastBends[activeIndex] = reverseBend;
                        Interval activeInterval = IKUtility.AngleInterval(prevBends[activeIndex], lastBends[activeIndex]);

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

            for (int j = 0; j < 2; j++)
            {
                for (int i = 0; i < angles[j].Count; i++)
                {
                    Debug.Log("i: " + j + " :" + angles[j][i]);
                }
            }
            for (int i = angles[1].Count - 1; i >= 0; i--)
            {
                this.angles[angles[0].Count + (angles[1].Count - 1 - i)] = angles[1][i];
            }
            if (shorterSide == IntervalSide.Min)
            {
                /*for (int i = 0; i < angles[1].Count; i++)
                {
                    this.angles[i] = angles[1][i];
                }
                for (int i = angles[0].Count - 1; i >= 0; i--)
                {
                    this.angles[angles[1].Count + (angles[0].Count - 1 - i)] = angles[0][i];
                }
                for (int i = 0; i < regions[1].Count; i++)
                {
                    this.regions[i] = regions[1][i];
                }
                for (int i = regions[0].Count - 1; i >= 0; i--)
                {
                    this.regions[regions[1].Count + (regions[0].Count - 1 - i)] = regions[0][i];
                }*/

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

        public float Get(IKSegment[] segments, float angle)
        {
            for (int i = 0; i < regions.Length; i++)
            {
                if (new Interval(angles[i], angles[i + 1]).Contains(angle))
                {
                    float angleSum = 0;
                    Vector2 center = Vector2.zero;
                    for (int j = startN; j < regions[i].MoveChainIndex; j++)
                    {
                        angleSum += segments[j].angleLimit[(int)regions[i].StartSide];
                        center += IKUtility.ToVector(angleSum, segments[j].length);
                    }
                    angleSum = 0;
                    Vector2 movingPart = Vector2.zero;
                    for (int j = regions[i].MoveChainIndex; j < partCount; j++)
                    {
                        angleSum += segments[j].angleLimit[(int)regions[i].ContinueSide];
                        movingPart += IKUtility.ToVector(angleSum, segments[j].length);
                    }

                    return IKUtility.IntersectionLineCircle(angle, center, movingPart.magnitude).p1.magnitude;
                }
            }
            return float.NaN;
        }

        public Circle GetCircle(IKSegment[] segments, int i)
        {
            Vector2 center = IKUtility.FollowChain(segments, new IndexRange(startN, regions[i].MoveChainIndex), regions[i].StartSide);
            float radius = IKUtility.FollowChain(segments, new IndexRange(regions[i].MoveChainIndex, partCount), regions[i].ContinueSide).magnitude;
            return new Circle(center, radius);
        }
        public Interval GetCircleInterval(IKSegment[] segments, int i)
        {
            Vector2 center = IKUtility.FollowChain(segments, new IndexRange(startN, regions[i].MoveChainIndex), regions[i].StartSide);
            Vector2 v1 = IKUtility.ToVector(angles[i], Get(segments, angles[i]));
            Vector2 v2 = IKUtility.ToVector(angles[i + 1], Get(segments, angles[i + 1]));
            return new Interval(IKUtility.ToAngle(v1 - center), IKUtility.ToAngle(v2 - center));
        }

        public int AngleCount { get { return angles.Length; } }
        public Interval ValidInterval { get { return new Interval(angles[0], angles[angles.Length - 1]); } }
    }
    [System.Serializable]
    public class MaxReachData
    {
        // from center angle interval those more outside
        private Interval[] angleSums;
        private int startN, endN;

        public MaxReachData(IKSegment[] segments, int startN, int endN)
        {
            this.startN = startN;
            this.endN = endN;

            angleSums = new Interval[endN - startN];
            for (int i = startN; i < endN; i++)
            {
                angleSums[i - startN] = IKUtility.AngleInterval(
                    IKUtility.FollowChain(segments, new IndexRangeK(startN, i + 1, endN), IntervalSide.Min, IntervalSide.Center),
                    IKUtility.FollowChain(segments, new IndexRangeK(startN, i + 1, endN), IntervalSide.Max, IntervalSide.Center));
            }
        }

        public float Get(IKSegment[] segments, float angle)
        {
            float[] lengths = new float[endN - startN];
            for (int i = startN; i < endN; i++)
            {
                lengths[i - startN] = segments[i].length;
            }
            for (int i = 0; i < angleSums.Length; i++)
            {
                //Debug.DrawAngleInterval(Vector2.zero, angleSums[i], 20, Color.magenta, 10);
            }
            return Get(segments, angle, endN, lengths);
        }
        private float Get(IKSegment[] segments, float angle, int k, float[] lengths)
        {
            // to prevent infinite loops (just in case though it shouldn't happen)
            if (k <= startN) { return 0; }
            // single chain part
            else if (k - startN == 1) { return lengths[0]; }

            if (angleSums[k - 2 - startN].Contains(angle))
            {
                float[] newLengths = new float[lengths.Length - 1];
                for (int i = 0; i < newLengths.Length; i++)
                {
                    newLengths[i] = lengths[i];
                }
                newLengths[newLengths.Length - 1] += lengths[newLengths.Length];
                return Get(segments, angle, k - 1, newLengths);
            }
            else if (new Interval(angleSums[k - 1 - startN].min, angleSums[k - 2 - startN].min).Contains(angle))
            {
                return IKUtility.IntersectionLineCircle(angle, IKUtility.FollowChain(segments, new IndexRange(startN, k - 1), IntervalSide.Min), lengths[k - 1 - startN]).p1.magnitude;
            }
            else if (new Interval(angleSums[k - 2 - startN].max, angleSums[k - 1 - startN].max).Contains(angle))
            {
                return IKUtility.IntersectionLineCircle(angle, IKUtility.FollowChain(segments, new IndexRange(startN, k - 1), IntervalSide.Max), lengths[k - 1 - startN]).p1.magnitude;
            }
            else
            {
                return float.NaN;
            }
        }

        public Circle GetCircle(IKSegment[] segments, int i)
        {
            int index = startN + Mathf.Abs(i); // from ... -2, -1, 0, 1, 2, ... to actual segment index
            Vector2 center = IKUtility.FollowChain(segments, new IndexRange(startN, index), (i < 0) ? IntervalSide.Min : IntervalSide.Max);
            float radius = 0;
            for (int k = index; k < endN; k++)
            {
                radius += segments[k].length;
            }
            return new Circle(center, radius);
        }
        public Interval GetCircleInterval(IKSegment[] segments, int i)
        {
            if (i == 0)
            {
                return angleSums[0];
            }
            else
            {
                IntervalSide side = (i < 0) ? IntervalSide.Min : IntervalSide.Max;
                int index = startN + Mathf.Abs(i); // from ... -2, -1, 0, 1, 2, ... to actual segment index
                Vector2 center = IKUtility.FollowChain(segments, new IndexRange(startN, index), side);
                Vector2 v1 = IKUtility.FollowChain(segments, new IndexRangeK(startN, index, endN), side, IntervalSide.Center);
                Vector2 v2 = IKUtility.FollowChain(segments, new IndexRangeK(startN, index + 1, endN), side, IntervalSide.Center);
                v1 = v1 - center;
                v2 = v2 - center;
                return IKUtility.AngleInterval(v1, v2);
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
    [System.Serializable]
    public class IKSegment
    {
        public float length;
        public Interval angleLimit;
        [Range(0, 1)] public float targetAngle = 0.5f;
        public bool animateTargetAngle = false;

        public IKSegment(float length, float angleMin, float angleMax)
        {
            this.length = length;
            this.angleLimit = new Interval(angleMin, angleMax);
        }
    }
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
    public struct Circle
    {
        public Vector2 center;
        public float radius;

        public Circle(Vector2 center, float radius)
        {
            this.center = center;
            this.radius = radius;
        }

        public override string ToString()
        {
            return $"({center}, {radius})";
        }
    }
    public struct CircleIntersection
    {
        public enum Variant { Miss, Touching, Intersect, Covering }

        public Vector2 p1, p2;

        public CircleIntersection(float x1, float y1, float x2, float y2) : this(new Vector2(x1, y1), new Vector2(x2, y2)) { }
        public CircleIntersection(Vector2 p1, Vector2 p2)
        {
            this.p1 = p1;
            this.p2 = p2;
        }

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
        public static CircleIntersection Miss { get { return new CircleIntersection(float.NaN, float.NaN, float.NaN, float.NaN); } }
    }
    public class IndexRange
    {
        public int startN, endN;

        public IndexRange(int startN, int partCount)
        {
            this.startN = startN;
            this.endN = partCount;
        }
    }
    public class IndexRangeK : IndexRange
    {
        public int k;

        public IndexRangeK(int startN, int k, int endN) : base(startN, endN)
        {
            this.k = k;
        }
    }

}