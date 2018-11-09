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
