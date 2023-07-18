using UnityEngine;
using UnityEditor;

namespace InverseKinematics2D
{
    [CustomEditor(typeof(IK2))]
    public class IK2Editor : Editor
    {
        public override void OnInspectorGUI()
        {
            base.OnInspectorGUI();
            if (GUILayout.Button("Refresh"))
            {
                ((IK2)target).Init();
            }
        }
    }

    [CustomEditor(typeof(IK2Segment)), CanEditMultipleObjects]
    public class IK2SegmentEditor : Editor
    {
        public override void OnInspectorGUI()
        {
            base.OnInspectorGUI();
            if (GUILayout.Button("Refresh"))
            {
                ((IK2Segment)target).GetComponentInParent<IK2>().Init();
            }
        }
    }
}