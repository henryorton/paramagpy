{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. automodule:: {{ fullname }}

   {% block functions %}
      {% if functions %}
         .. rubric:: Functions

         .. autosummary::
            :toctree:
            {% for item in functions %}
               {{ item }}
            {%- endfor %}
      {% endif %}
   {% endblock %}

   {% block classes %}
      {% if classes %}
         .. rubric:: Classes

         .. autosummary::
            :toctree:
            {% for item in classes %}
               {{ item }}
            {%- endfor %}
      {% endif %}
   {% endblock %}
