using UnityEditor;
using UnityEngine;

namespace InverseKinematics1
{
    [CustomPropertyDrawer(typeof(Constraint))]
    public class ConstraintDrawer : PropertyDrawer
    {
        public override void OnGUI(Rect position, SerializedProperty property, GUIContent label)
        {
            EditorGUI.BeginProperty(position, label, property);

            float temp = EditorGUIUtility.labelWidth;
            EditorGUIUtility.labelWidth = 25.0f;

            EditorGUI.LabelField(new Rect(position.x, position.y, position.width / 3, position.height), label);
            EditorGUI.PropertyField(new Rect(position.x + position.width / 3, position.y, position.width / 3, position.height), property.FindPropertyRelative("min"));
            EditorGUI.PropertyField(new Rect(position.x + position.width * 2 / 3, position.y, position.width / 3, position.height), property.FindPropertyRelative("max"));

            EditorGUIUtility.labelWidth = temp;

            EditorGUI.EndProperty();
        }
    }

}