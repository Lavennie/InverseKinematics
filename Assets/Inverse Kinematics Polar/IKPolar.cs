using UnityEngine;
using System.Collections.Generic;
#if UNITY_EDITOR
using UnityEditor;
#endif

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
        Vector2 p2 = FollowChain(segments, new IndexRange(fromToSplit.k, fromToSplit.partCount), limitIndex2);
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
        for (int i = fromTo.startN; i < fromTo.partCount; i++)
        {
            if (limitIndex == IntervalSide.Min || limitIndex == IntervalSide.Max)
            {
                endAngle += segments[i].angleLimit[(int)limitIndex];
            }
            p += ToVector(endAngle, segments[i].length);
        }
        return p;
    }

    public static void DebugFollowChain(IKSegment[] segments, IndexRangeK fromToSplit, IntervalSide limitIndex1, IntervalSide limitIndex2, float startAngle)
    {
        Vector2 p1 = FollowChain(segments, new IndexRange(fromToSplit.startN, fromToSplit.k), limitIndex1, out float endAngle);
        DebugFollowChain(segments, new IndexRange(fromToSplit.startN, fromToSplit.k), limitIndex1, startAngle, 
            FollowChain(segments, new IndexRange(0, fromToSplit.startN), limitIndex1), out float endAngle1);
        DebugFollowChain(segments, new IndexRange(fromToSplit.k, fromToSplit.partCount), limitIndex2, endAngle1, p1);
    }
    public static void DebugFollowChain(IKSegment[] segments, IndexRange fromTo, IntervalSide limitIndex, float startAngle, Vector2 origin)
    {
        DebugFollowChain(segments, fromTo, limitIndex, startAngle, origin, out float endAngle);
    }
    public static void DebugFollowChain(IKSegment[] segments, IndexRange fromTo, IntervalSide limitIndex, float startAngle, Vector2 origin, out float endAngle)
    {
        endAngle = startAngle;
        Vector2 p = origin;
        for (int i = fromTo.startN; i < fromTo.partCount; i++)
        {
            if (limitIndex == IntervalSide.Min || limitIndex == IntervalSide.Max)
            {
                endAngle += segments[i].angleLimit[(int)limitIndex];
            }
            Vector2 nextChain = ToVector(endAngle, segments[i].length);
            Debug.DrawLine(p, p + nextChain, Color.yellow, 10);
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
        float x1 = center.x + Mathf.Sqrt(radius * radius - (y1 - center.y) * (y1 - center.y));
        float x2 = center.x - Mathf.Sqrt(radius * radius - (y2 - center.y) * (y2 - center.y));

        return new CircleIntersection(x1, y1, x2, y2);
    }
}

public class IKPolar : MonoBehaviour
{
    public bool debugChainBend = true;
    public bool debugTargetLine = true;
    public Transform target;

    public IKSegment[] segments;
    private Reach[] dataFront;
    private Reach[] dataBack;

    private Vector2 mousePos = Vector2.zero;
    private Vector2 reach = Vector2.zero;

    private void Awake()
    {
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
            DebugChainBend(Vector2.zero, 0, 0, Color.green);
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

        Debug.DrawReach(dataFront[0], segments, Color.green);
        //Align2(Vector3.zero, mousePos, -45, segments.Length);
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

    private void DebugChainBend(Vector2 pos, float angle, int startN, Color color)
    {
        DebugChainBend(pos, Quaternion.Euler(0, 0, angle), startN, color);
    }
    private void DebugChainBend(Vector2 pos, Quaternion rot, int startN, Color color)
    {
        // draw min/max bend
        Vector2 pMin = pos;
        Vector2 pMax = pos;
        Interval angleSum = new Interval(0, 0);
        for (int i = 0; i < segments.Length; i++)
        {
            angleSum.min += segments[i].angleLimit[0];
            angleSum.max += segments[i].angleLimit[1];
            Vector2 dirMin = rot * new Vector2(Mathf.Cos(angleSum.min * Mathf.Deg2Rad), Mathf.Sin(angleSum.min * Mathf.Deg2Rad)) * segments[i].length;
            Vector2 dirMax = rot * new Vector2(Mathf.Cos(angleSum.max * Mathf.Deg2Rad), Mathf.Sin(angleSum.max * Mathf.Deg2Rad)) * segments[i].length;
            Debug.DrawRay(pMin, dirMin, color, 0.01f);
            Debug.DrawRay(pMax, dirMax, color, 0.01f);
            pMin += dirMin;
            pMax += dirMax;
        }
    }

    private void Precalculate()
    {
        dataFront = new Reach[segments.Length];
        int i = 0;
        for (int k = segments.Length; k > segments.Length - 1; k--)
        {
            dataFront[i] = new Reach(segments, 0, k);
            i++;
        }

        /*dataBack = new Reach[segments.Length];
        for (int j = 0; j < segments.Length; j++)
        {
            dataBack[j] = new Reach(segments, j, segments.Length);
        }*/
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




    public void Align2(Vector3 origin, Vector2 target, float angle, int n)
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

        Debug.Log(n);
        int segmentI = segments.Length - n;
        Debug.DrawCircle(origin, Vector3.Distance(origin, target), 64, Color.black);
        DebugChainBend(origin, segments[segmentI].angleLimit.min, segmentI, Color.red);
        DebugChainBend(origin, segments[segmentI].angleLimit.max, segmentI, Color.magenta);

        Debug.DrawReach(dataBack[segmentI], segments, Color.red, segments[segmentI].angleLimit.min);
        Debug.DrawReach(dataBack[segmentI], segments, Color.magenta, segments[segmentI].angleLimit.max);
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

    public MinReachData(IKSegment[] segments, int startN, int partCount)
    {
        this.startN = startN;
        this.partCount = partCount;

        Debug.Log("_ " + startN + "  " + partCount);

        if (partCount - startN <= 0)
        {
            return;
        }
        // only one chain part
        if (partCount - startN == 1)
        {
            this.angles = new float[2] { segments[startN].angleLimit.min, segments[startN].angleLimit.max };
            this.regions = new MinReachRegion[1] { new MinReachRegion(0, 1, startN) };
            return;
        }
        // [0]: shorter bend side
        // [1]: longer bend side
        int shorterIndex = 0;
        Vector3[] bendPoints = new Vector3[2]
        {
            IKUtility.FollowChain(segments, new IndexRange(startN, partCount), IntervalSide.Min),
            IKUtility.FollowChain(segments, new IndexRange(startN, partCount), IntervalSide.Max)
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
            shorterIndex = 1;
        }

        // chain with two parts
        if (partCount - startN == 2)
        {
            this.angles = new float[3]
            {
                IKUtility.ToAngle(IKUtility.FollowChain(segments, new IndexRange(startN, partCount), (IntervalSide)shorterIndex)),
                IKUtility.ToAngle(IKUtility.FollowChain(segments, new IndexRangeK(startN, startN + 1, partCount), (IntervalSide)((shorterIndex + 1) % 2), (IntervalSide)shorterIndex)),
                IKUtility.ToAngle(IKUtility.FollowChain(segments, new IndexRange(startN, partCount), (IntervalSide)((shorterIndex + 1) % 2)))
            };
            this.regions = new MinReachRegion[2] { new MinReachRegion((shorterIndex + 1) % 2, shorterIndex, startN), new MinReachRegion((shorterIndex + 1) % 2, shorterIndex, startN + 1) };
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
        while (true)
        {
            bool exit = false;
            // radius of above interval at end point that is more bent
            float activeLength;
            float inactiveLength;
            do
            {
                IntervalSide activeSide = (activeIndex == 0) ? (IntervalSide)shorterIndex : (IntervalSide)((shorterIndex + 1) % 2);
                IntervalSide inactiveSide = (activeSide == IntervalSide.Min) ? IntervalSide.Max : IntervalSide.Min;

                // calculate the bend position (Forward Kinematics)
                IKUtility.DebugFollowChain(segments, new IndexRangeK(startN, k[activeIndex], partCount), inactiveSide, activeSide, 0);
                Vector2 reverseBend = IKUtility.FollowChain(segments, new IndexRangeK(startN, k[activeIndex], partCount), inactiveSide, activeSide);
                Vector2 debugCenter = IKUtility.FollowChain(segments, new IndexRange(0, k[activeIndex] - 1), inactiveSide);
                float debugRadius = IKUtility.FollowChain(segments, new IndexRange(k[activeIndex] - 1, partCount), activeSide).magnitude;
                Debug.DrawCircle(debugCenter, debugRadius, 64, Color.cyan, 10);
                //Debug.Log(activeSide + "  " + k[activeIndex]);
                // update bend history
                prevBends[activeIndex] = lastBends[activeIndex];
                lastBends[activeIndex] = reverseBend;
                // arc angle interval from prevBend to lastBend
                Interval activeInterval = IKUtility.AngleInterval(prevBends[activeIndex], lastBends[activeIndex]);
                Interval inactiveInterval = IKUtility.AngleInterval(prevBends[inactiveIndex], lastBends[inactiveIndex]);

                //Debug.DrawAngleInterval(debugCenter, activeInterval, debugRadius, Color.white, 10);
                //Debug.DrawAngleInterval(debugCenter, inactiveInterval, debugRadius, Color.black, 10);
                if (Interval.AreOverlapping(activeInterval, inactiveInterval))
                {
                    Vector2 c1 = IKUtility.FollowChain(segments, new IndexRange(startN, k[activeIndex]), inactiveSide);
                    Vector2 c2 = IKUtility.FollowChain(segments, new IndexRange(startN, k[inactiveIndex] - 1), activeSide);
                    Vector2 r1 = IKUtility.FollowChain(segments, new IndexRange(k[activeIndex], partCount), activeSide);
                    Vector2 r2 = IKUtility.FollowChain(segments, new IndexRange(k[inactiveIndex] - 1, partCount), inactiveSide);
                    //Debug.DrawCircle(c1, r1.magnitude, 128, Color.white, 10);
                    //Debug.DrawCircle(c2, r2.magnitude, 128, Color.black, 10);
                    CircleIntersection intersection = IKUtility.TwoCircleIntersection(c1, c2, r1.magnitude, r2.magnitude);

                    if (IKUtility.AngleInterval(prevBends[inactiveIndex], lastBends[activeIndex]).Contains(IKUtility.ToAngle(intersection.p1)))
                    {
                        angles[activeIndex].Insert(insertI[activeIndex], IKUtility.ToAngle(lastBends[activeIndex]));
                        regions[activeIndex].Insert(insertI[activeIndex], new MinReachRegion((int)inactiveSide, (int)activeSide, k[activeIndex]));
                        angles[activeIndex].Insert(insertI[activeIndex], IKUtility.ToAngle(intersection.p1));
                        regions[activeIndex].Insert(insertI[activeIndex], new MinReachRegion((int)activeSide, (int)inactiveSide, k[inactiveIndex] - 1));
                        insertI[activeIndex]++;
                    }
                    else
                    {
                        // TODO:
                    }
                    exit = true;
                    break;
                }
                // continue bending the same side
                else
                {
                    Debug.Log("insert " + activeIndex);
                    angles[activeIndex].Insert(insertI[activeIndex], activeInterval[inactiveIndex]);
                    regions[activeIndex].Insert(insertI[activeIndex], new MinReachRegion((int)inactiveSide, (int)activeSide, k[activeIndex] - 1));
                    insertI[activeIndex]++;
                    k[activeIndex]++;
                }

                activeLength = lastBends[activeIndex].magnitude;
                inactiveLength = lastBends[inactiveIndex].magnitude;

                // infinite loop protection while still writing this function
                if (k[0] > 100 || k[1] > 100)
                {
                    exit = true;
                    break;
                }
            } while (activeLength < inactiveLength); // not if at end activeLength < inactiveLength, but if at any time it crossed over inactiveLength
            activeIndex = inactiveIndex; 
            inactiveIndex = (inactiveIndex + 1) % 2;

            if (exit)
            {
                break;
            }
        }

        // add angle to full bend on shorter side
        angles[shorterIndex].Insert(0, IKUtility.ToAngle(bendPoints[shorterIndex]));
        // add angle to full bend on longer side
        angles[(shorterIndex + 1) % 2].Add(IKUtility.ToAngle(bendPoints[(shorterIndex + 1) % 2]));
        regions[(shorterIndex + 1) % 2].Add(new MinReachRegion(shorterIndex, (shorterIndex + 1) % 2, startN));

        for (int i = 0; i < regions[0].Count; i++)
        {
            Debug.Log("shorter: " + regions[0][i]);
        }
        for (int i = 0; i < regions[1].Count; i++)
        {
            Debug.Log("longer: " + regions[1][i]);
        }

        // save calculated values
        this.angles = new float[angles[0].Count + angles[1].Count];
        for (int i = 0; i < angles[0].Count; i++)
        {
            this.angles[i] = angles[0][i];
            Debug.DrawAngle(Vector3.zero, this.angles[i], 10, Color.black, 10);
        }
        for (int i = 0; i < angles[1].Count; i++)
        {
            this.angles[angles[0].Count + i] = angles[1][i];
            Debug.DrawAngle(Vector3.zero, this.angles[angles[0].Count + i], 10, Color.black, 10);
        }
        this.regions = new MinReachRegion[regions[0].Count + regions[1].Count];
        for (int i = 0; i < regions[0].Count; i++)
        {
            this.regions[i] = regions[0][i];
        }
        for (int i = 0; i < regions[1].Count; i++)
        {
            this.regions[regions[0].Count + i] = regions[1][i];
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
                for (int j = 0; j < regions[i].MoveChainIndex; j++)
                {
                    angleSum += segments[j].angleLimit[regions[i].StartSide];
                    center += IKUtility.ToVector(angleSum, segments[j].length);
                }
                angleSum = 0;
                Vector2 movingPart = Vector2.zero;
                for (int j = regions[i].MoveChainIndex; j < partCount; j++)
                {
                    angleSum += segments[j].angleLimit[regions[i].ContinueSide];
                    movingPart += IKUtility.ToVector(angleSum, segments[j].length);
                }

                //Debug.DrawCircle(new Circle(center, movingPart.magnitude), 128, Color.yellow);

                return IKUtility.IntersectionLineCircle(angle, center, movingPart.magnitude).p1.magnitude;
            }
        }
        return float.NaN;
    }

    public Circle GetCircle(IKSegment[] segments, int i)
    {
        Vector2 center = IKUtility.FollowChain(segments, new IndexRange(startN, regions[i].MoveChainIndex), (IntervalSide)regions[i].StartSide);
        float radius = IKUtility.FollowChain(segments, new IndexRange(regions[i].MoveChainIndex, partCount), (IntervalSide)regions[i].ContinueSide).magnitude;
        return new Circle(center, radius);
    }
    public Interval GetCircleInterval(IKSegment[] segments, int i)
    {
        Vector2 center = IKUtility.FollowChain(segments, new IndexRange(startN, regions[i].MoveChainIndex), (IntervalSide)regions[i].StartSide);
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
        float[] lengths = new float[endN];
        for (int i = startN; i < endN; i++)
        {
            lengths[i] = segments[i].length;
        }
        return Get(segments, angle, endN, lengths);
    }
    private float Get(IKSegment[] segments, float angle, int k, float[] lengths)
    {
        // to prevent infinite loops (just in case though it shouldn't happen)
        if (k <= startN) { return 0; }
        // single chain part
        else if (k - startN == 1) { return lengths[startN]; }

        if (angleSums[k - 2].Contains(angle))
        {
            float[] newLengths = new float[lengths.Length - 1];
            for (int i = 0; i < newLengths.Length; i++)
            {
                newLengths[i] = lengths[i];
            }
            newLengths[newLengths.Length - 1] += lengths[newLengths.Length];
            return Get(segments, angle, k - 1, newLengths);
        }
        else if (new Interval(angleSums[k - 1].min, angleSums[k - 2].min).Contains(angle))
        {
            return IKUtility.IntersectionLineCircle(angle, IKUtility.FollowChain(segments, new IndexRange(startN, k - 1), IntervalSide.Min), lengths[k - 1]).p1.magnitude;
        }
        else if (new Interval(angleSums[k - 2].max, angleSums[k - 1].max).Contains(angle))
        {
            return IKUtility.IntersectionLineCircle(angle, IKUtility.FollowChain(segments, new IndexRange(startN, k - 1), IntervalSide.Max), lengths[k - 1]).p1.magnitude;
        }
        else
        {
            return float.NaN;
        }
    }

    public Circle GetCircle(IKSegment[] segments, int i)
    {
        int index = Mathf.Abs(i); // from ... -2, -1, 0, 1, 2, ... to actual segment index
        Vector2 center = IKUtility.FollowChain(segments, new IndexRange(startN, index), (i < 0) ? IntervalSide.Min : IntervalSide.Max);
        float radius = 0;
        for (int k = index; k < endN; k++)
        {
            radius += segments[startN + k].length;
        }
        return new Circle(center, radius);
    }
    public Interval GetCircleInterval(IKSegment[] segments, int i)
    {
        if (i == 0)
        {
            return angleSums[startN];
        }
        else
        {
            IntervalSide side = (i < 0) ? IntervalSide.Min : IntervalSide.Max;
            int index = Mathf.Abs(i); // from ... -2, -1, 0, 1, 2, ... to actual segment index
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
    private int startSide;
    private int continueSide;
    private int moveChainIndex;

    public MinReachRegion(int startSide, int continueSide, int moveChainIndex)
    {
        this.startSide = startSide;
        this.continueSide = continueSide;
        this.moveChainIndex = moveChainIndex;
    }

    public override string ToString()
    {
        return $"(Start side = {startSide}, continue side = {continueSide}, move chain index = {moveChainIndex})";
    }

    public int StartSide { get { return startSide; } }
    public int ContinueSide { get { return continueSide; } }
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
    public int startN, partCount;

    public IndexRange(int startN, int partCount)
    {
        this.startN = startN;
        this.partCount = partCount;
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