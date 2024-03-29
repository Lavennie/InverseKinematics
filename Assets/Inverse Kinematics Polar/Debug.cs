using UnityEngine;
using InverseKinematicsPolarTest;

public class Debug : UnityEngine.Debug
{
    public static void DrawReach(Reach reach, IKSegment[] segments, Color color, float rotAngle = 0)
    {
        DrawReach(Vector3.zero, reach, segments, color, 0, rotAngle);
    }
    public static void DrawReach(Vector3 origin, Reach reach, IKSegment[] segments, Color color, float rotAngle = 0)
    {
        DrawReach(origin, reach, segments, color, 0, rotAngle);
    }
    public static void DrawReach(Vector3 origin, Reach reach, IKSegment[] segments, Color color, float duration, float rotAngle = 0)
    {
        foreach (var circleInterval in reach.GetMaxCircles(segments))
        {
            Circle circle = new Circle(origin + Quaternion.Euler(0, 0, rotAngle) * circleInterval.Key.center, circleInterval.Key.radius);
            DrawArc(circle, circleInterval.Value, 64, color, duration, rotAngle);
        }
        color *= 0.7f;
        color.a = 1.0f;
        foreach (var circleInterval in reach.GetMinCircles(segments))
        {
            Circle circle = new Circle(origin + Quaternion.Euler(0, 0, rotAngle) * circleInterval.Key.center, circleInterval.Key.radius);
            DrawArc(circle, circleInterval.Value, 64, color, duration, rotAngle);
        }
    }

    public static void DrawCircle(Circle circle, int segments, Color color)
    {
        DrawCircle(circle.center, circle.radius, segments, color);
    }
    public static void DrawCircle(Circle circle, int segments, Color color, float duration)
    {
        DrawCircle(circle.center, circle.radius, segments, color, duration);
    }
    public static void DrawCircle(Vector3 position, float radius, int segments, Color color)
    {
        DrawCircle(position, radius, segments, color, 0);
    }
    public static void DrawCircle(Vector3 position, float radius, int segments, Color color, float duration)
    {
        if (radius <= 0.0f || segments <= 0) { return; }

        float angleStep = (360.0f / segments) * Mathf.Deg2Rad;

        Vector3 lineStart;
        Vector3 lineEnd;
        for (int i = 0; i < segments; i++)
        {
            lineStart = position + new Vector3(Mathf.Cos(angleStep * i), Mathf.Sin(angleStep * i)) * radius;
            lineEnd = position + new Vector3(Mathf.Cos(angleStep * (i + 1)), Mathf.Sin(angleStep * (i + 1))) * radius;

            DrawLine(lineStart, lineEnd, color, duration);
        }

        DrawPoint(position, color, duration, radius / 10);
    }

    public static void DrawArc(Circle circle, Interval angle, int segments, Color color, float rotAngle = 0)
    {
        DrawArc(circle.center, circle.radius, angle, segments, color, 0, rotAngle);
    }
    public static void DrawArc(Circle circle, Interval angle, int segments, Color color, float duration, float rotAngle = 0)
    {
        DrawArc(circle.center, circle.radius, angle, segments, color, duration, rotAngle);
    }
    public static void DrawArc(Vector3 center, float radius, Interval angle, int segments, Color color, float rotAngle = 0)
    {
        DrawArc(center, radius, angle, segments, color, 0, rotAngle);
    }
    public static void DrawArc(Vector3 center, float radius, Interval angle, int segments, Color color, float duration, float rotAngle = 0)
    {
        if (radius <= 0.0f || segments <= 0) { return; }

        float angleStep = ((angle.max - angle.min) / segments);
        angleStep *= Mathf.Deg2Rad;
        angle.min *= Mathf.Deg2Rad;
        angle.max *= Mathf.Deg2Rad;
        Vector3 lineStart;
        Vector3 lineEnd;
        for (int i = 0; i < segments; i++)
        {
            lineStart = center + Quaternion.Euler(0, 0, rotAngle) * new Vector3(Mathf.Cos(angle.min + angleStep * i), Mathf.Sin(angle.min + angleStep * i)) * radius;
            lineEnd = center + Quaternion.Euler(0, 0, rotAngle) * new Vector3(Mathf.Cos(angle.min + angleStep * (i + 1)), Mathf.Sin(angle.min + angleStep * (i + 1))) * radius;

            DrawLine(lineStart, lineEnd, color, duration);
        }

        DrawPoint(center, color, duration, radius / 10);
    }

    public static void DrawAngle(Vector3 position, float angle, float radius, Color color)
    {
        DrawRay(position, new Vector2(Mathf.Cos(angle * Mathf.Deg2Rad), Mathf.Sin(angle * Mathf.Deg2Rad)) * radius, color);
    }
    public static void DrawAngle(Vector3 position, float angle, float radius, Color color, float duration)
    {
        DrawRay(position, new Vector2(Mathf.Cos(angle * Mathf.Deg2Rad), Mathf.Sin(angle * Mathf.Deg2Rad)) * radius, color, duration);
    }

    public static void DrawAngleInterval(Vector3 position, Interval angleInterval, float radius, Color color)
    {
        DrawAngle(position, angleInterval.min, radius, color);
        DrawAngle(position, angleInterval.max, radius, color);
        DrawArc(position, radius / 6, angleInterval, 64, color);
    }
    public static void DrawAngleInterval(Vector3 position, Interval angleInterval, float radius, Color color, float duration)
    {
        DrawAngle(position, angleInterval.min, radius, color, duration);
        DrawAngle(position, angleInterval.max, radius, color, duration);
        DrawArc(position, radius / 6, angleInterval, 64, color, duration);
    }

    public static void DrawPoint(Vector3 position, Color color, float size = 0.1f)
    {
        DrawLine(position - new Vector3(size, 0, 0), position + new Vector3(size, 0, 0), color);
        DrawLine(position - new Vector3(0, size, 0), position + new Vector3(0, size, 0), color);
    }
    public static void DrawPoint(Vector3 position, Color color, float duration, float size = 0.1f)
    {
        DrawLine(position - new Vector3(size, 0, 0), position + new Vector3(size, 0, 0), color, duration);
        DrawLine(position - new Vector3(0, size, 0), position + new Vector3(0, size, 0), color, duration);
    }

    public static void DrawPoints(CircleIntersection intersection, Color color, float size = 0.1f)
    {
        switch (intersection.Type)
        {
            case CircleIntersection.Variant.Touching:
                DrawPoint(intersection.p1, color, size);
                break;
            case CircleIntersection.Variant.Intersect:
                DrawPoint(intersection.p1, color, size);
                DrawPoint(intersection.p2, color, size);
                break;
        }
    }
    public static void DrawPoints(CircleIntersection intersection, Color color, float duration, float size = 0.1f)
    {
        switch (intersection.Type)
        {
            case CircleIntersection.Variant.Touching:
                DrawPoint(intersection.p1, color, duration, size);
                break;
            case CircleIntersection.Variant.Intersect:
                DrawPoint(intersection.p1, color, duration, size);
                DrawPoint(intersection.p2, color, duration, size);
                break;
        }
    }
}
