define([
    'app'
], function (App) {

    App.View.List = App.Backbone.View.extend({
        tagName: 'div',
        className: 'thumb-list',

        $emptyMessage: null,

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

        showEmptyMessage: function(msg) {
            msg = msg || 'No result';

            if (this.subViews && this.subViews.length) {
                this.hideEmptyMessage();
            } else {
                this.$emptyMessage.html(msg).show();
            }

            return this;
        },

        hideEmptyMessage: function() {
            this.$emptyMessage.hide();

            return this;
        },

        render: function () {
            this.$emptyMessage = $('<p></p>')
                .addClass('list-empty-message')
                .appendTo(this.$el)
                .hide();

            return this;
        }
    });

    App.View.List.create = function(options, createOptions) {
        options = options || {};

        var view = new App.View.List(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.List;
});
