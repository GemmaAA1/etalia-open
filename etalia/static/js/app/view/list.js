define([
    'app'
], function (App) {

    App.View.List = App.Backbone.View.extend({
        tagName: 'div',
        className: 'thumb-list',

        addThumbView: function(view) {
            this.$el.append(view.render().$el);
            this.pushSubView(view);
        },

        removeThumbById: function(id) {
            if (this.subViews) {
                var filter = function (view) {
                        return view.id === id;
                    },
                    view = App._.find(this.subViews, filter);
                if (view) {
                    view.remove();
                    this.subViews = App._.reject(this.subViews, filter);
                }
            }
        },

        clear: function() {
            if (this.subViews) {
                App._.each(this.subViews, function(view) {
                   view.remove();
                });
            }
            this.subViews = null;
        },

        render: function () {
            //this.$el.html(this.template({}));

            return this;
        }
    });

    App.View.List.create = function(options, createOptions) {
        options = options || {};
        /*if (!options.tabs) {
            throw 'options.tabs is expected to be a hash of tabs';
        }*/

        var view = new App.View.List(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.List;
});
